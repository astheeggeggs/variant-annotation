#!/usr/bin/env python
"""
Helpers that sit between VEP and brava_create_annot.py.

  vep-tables      Rebuild the original pipeline's processed VEP tables (same 10 columns), filling
                  REVEL_SCORE / CADD_PHRED from the slim dbNSFP table exactly as the VEP dbNSFP
                  plugin would have, and list the variants that still need SpliceAI / CADD.
  cadd-prescored  Look up CADD v1.6 prescored indel scores; write the hits in CADD output format
                  and the misses as a VCF ready for CADD.sh.
"""

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "SAIGE_annotations" / "scripts"))
from brava_create_annot import (CADD_CUTOFF, MISSENSE_CSQS, REVEL_CUTOFF,  # noqa: E402
                                first_non_missing_revel, has_any_csq)

RAW_COLS = ["SNP_ID", "GENE", "LOF", "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL", "BIOTYPE", "AMINO_ACIDS"]
OUT_COLS = ["SNP_ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ", "TRANSCRIPT", "MANE_SELECT",
            "CANONICAL", "BIOTYPE"]
# Consequences the VEP 105 dbNSFP plugin annotates by default (%INCLUDE_SO in dbNSFP.pm)
DBNSFP_CSQS = ["missense_variant", "stop_lost", "stop_gained", "start_lost"]
# Latin-1 round-trips every byte, so tables are reproduced byte-for-byte whatever their encoding
ENC = "latin-1"


def log(msg):
    print(f"[brava_prep] {msg}", file=sys.stderr, flush=True)


def tabix_query(path, regions, columns):
    """Fetch rows overlapping `regions` (DataFrame chrom/pos) from a tabix-indexed file or URL."""
    if regions.empty:
        return pd.DataFrame(columns=columns)
    reg_file = Path(str(regions.attrs["tmp"]))
    regions.drop_duplicates().to_csv(reg_file, sep="\t", header=False, index=False)
    proc = subprocess.run(["tabix", "-R", str(reg_file), path], check=True, capture_output=True)
    reg_file.unlink()
    lines = [l.split("\t") for l in proc.stdout.decode(ENC).splitlines() if l and not l.startswith("#")]
    return pd.DataFrame([l[:len(columns)] for l in lines], columns=columns)


def split_snp_id(snp_id):
    parts = snp_id.str.split(":", n=3, expand=True)
    parts.columns = ["chrom", "pos", "ref", "alt"]
    return parts


def kept_rows(df):
    """protein_coding rows on MANE Select, or canonical where MANE is absent (as brava_create_annot)."""
    mane_missing = df["MANE_SELECT"] == "."
    return (df["BIOTYPE"] == "protein_coding") & (~mane_missing | (df["CANONICAL"] == "YES"))


def to_float(s):
    return pd.to_numeric(s.where(s != "."), errors="coerce")


def fill_dbnsfp(df, dbnsfp, tmpdir):
    """Reproduce the VEP 105 dbNSFP plugin (pep_match=1) for REVEL_score and CADD_phred."""
    v = split_snp_id(df["SNP_ID"])
    snv = v["ref"].str.fullmatch("[ACGT]") & v["alt"].str.fullmatch("[ACGT]")
    cand = snv & has_any_csq(df["CSQ"], DBNSFP_CSQS) & (df["AMINO_ACIDS"] != ".")
    df["REVEL_SCORE"] = "."
    df["CADD_PHRED"] = "."
    if not cand.any():
        return df
    # dbNSFP uses bare chromosome names and 'M' for the mitochondrion (as the plugin maps them)
    chrom = v.loc[cand, "chrom"].str.replace("^chr", "", regex=True).replace({"MT": "M"})
    regions = pd.DataFrame({"chrom": chrom, "beg": v.loc[cand, "pos"], "end": v.loc[cand, "pos"]})
    regions.attrs["tmp"] = Path(tmpdir) / "dbnsfp_regions.tsv"
    hits = tabix_query(dbnsfp, regions, ["chrom", "pos", "ref", "alt", "aa", "REVEL_SCORE", "CADD_PHRED"])
    log(f"dbNSFP: {cand.sum()} candidate rows, {len(hits)} table rows at those positions")
    # The table is already reduced to the first dbNSFP row per (chrom, pos, alt, aa), which is the
    # row the plugin's `last` picks; the plugin ignores REF.
    hits = hits.drop_duplicates(["chrom", "pos", "alt", "aa"], keep="first")
    key = pd.DataFrame({"chrom": chrom, "pos": v.loc[cand, "pos"], "alt": v.loc[cand, "alt"],
                        "aa": df.loc[cand, "AMINO_ACIDS"]})
    merged = key.reset_index().merge(hits.drop(columns="ref"), on=["chrom", "pos", "alt", "aa"],
                                     how="left").set_index("index")
    df.loc[merged.index, "REVEL_SCORE"] = merged["REVEL_SCORE"].fillna(".")
    df.loc[merged.index, "CADD_PHRED"] = merged["CADD_PHRED"].fillna(".")
    return df


def write_table(df, path):
    with open(path, "w", encoding=ENC, newline="") as f:
        f.write(" ".join(OUT_COLS) + "\n")
        df[OUT_COLS].to_csv(f, sep=" ", header=False, index=False)


def vep_tables(args):
    df = pd.read_csv(args.raw, sep=" ", header=None, names=RAW_COLS, dtype=str,
                     keep_default_na=False, encoding=ENC)
    log(f"{len(df)} VEP rows")
    df = fill_dbnsfp(df, args.dbnsfp, Path(args.out).parent)

    write_table(df, f"{args.out}.vep_processed.txt")

    with gzip.open(args.popmax, "rt") as f:
        common = set(line.strip() for line in f)
    df = df[~df["SNP_ID"].isin(common)]
    write_table(df, f"{args.out}.vep.gnomad_popmax_0.01_processed.txt")
    log(f"{len(df)} rows after gnomAD popmax > 0.01 filter")

    # Which rows can SpliceAI / CADD still change? Everything else is decided without them.
    kept = df[kept_rows(df)].copy()
    lof_hc = kept["LOF"] == "HC"
    missense = has_any_csq(kept["CSQ"], MISSENSE_CSQS)
    revel = kept["REVEL_SCORE"].map(lambda x: first_non_missing_revel(None if x == "." else x)).astype(float)
    cadd = to_float(kept["CADD_PHRED"])
    decided = lof_hc | (missense & ((revel >= REVEL_CUTOFF) | (cadd >= CADD_CUTOFF)))

    if args.all_spliceai:
        need_splice = pd.Series(kept["SNP_ID"].unique())
    else:
        need_splice = pd.Series(kept.loc[~decided, "SNP_ID"].unique())
    need_splice.to_csv(f"{args.out}.spliceai_ids.txt", header=False, index=False)
    log(f"SpliceAI needed for {len(need_splice)} variants "
        f"({kept['SNP_ID'].nunique()} with a MANE/canonical protein_coding row)")

    v = split_snp_id(kept["SNP_ID"])
    indel = v["ref"].str.len() != v["alt"].str.len()
    need_cadd = indel & missense & ~lof_hc & cadd.isna()
    sites = v[need_cadd].drop_duplicates()
    sites.to_csv(f"{args.out}.cadd_indel_sites.tsv", sep="\t", header=False, index=False)
    log(f"CADD needed for {len(sites)} indels (missense-class, not HC pLoF)")


def cadd_prescored(args):
    sites = pd.read_csv(args.sites, sep="\t", header=None, names=["chrom", "pos", "ref", "alt"], dtype=str)
    sites["Chrom"] = sites["chrom"].str.replace("^chr", "", regex=True)
    header = "## CADD GRCh38-v1.6 (c) University of Washington, Hudson-Alpha Institute for Biotechnology " \
             "and Berlin Institute of Health 2013-2020. All rights reserved.\n" \
             "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n"
    hits = pd.DataFrame(columns=["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"])
    if not sites.empty and args.prescored:
        regions = pd.DataFrame({"chrom": sites["Chrom"], "beg": sites["pos"], "end": sites["pos"]})
        regions.attrs["tmp"] = Path(args.out_tsv).parent / "cadd_regions.tsv"
        hits = tabix_query(args.prescored, regions, ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"])
        hits = hits.merge(sites.rename(columns={"pos": "Pos", "ref": "Ref", "alt": "Alt"})[
            ["Chrom", "Pos", "Ref", "Alt"]], on=["Chrom", "Pos", "Ref", "Alt"])
    # user-supplied CADD output (e.g. from CADD.sh on the unscored indels) takes precedence
    extra = [pd.read_csv(f, sep="\t", comment="#", header=None, dtype=str,
                         names=["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"]) for f in args.extra]
    if extra and not sites.empty:
        extra = pd.concat(extra)
        extra["Chrom"] = extra["Chrom"].str.replace("^chr", "", regex=True)
        extra = extra.merge(sites.rename(columns={"pos": "Pos", "ref": "Ref", "alt": "Alt"})[
            ["Chrom", "Pos", "Ref", "Alt"]], on=["Chrom", "Pos", "Ref", "Alt"])
        hits = pd.concat([extra, hits])
    hits = hits.drop_duplicates(["Chrom", "Pos", "Ref", "Alt"], keep="first")
    with gzip.open(args.out_tsv, "wt") as f:
        f.write(header)
        hits.to_csv(f, sep="\t", header=False, index=False)

    scored = set(zip(hits["Chrom"], hits["Pos"], hits["Ref"], hits["Alt"]))
    missing = sites[[k not in scored for k in zip(sites["Chrom"], sites["pos"], sites["ref"], sites["alt"])]]
    with open(args.out_unscored, "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in missing.itertuples():
            f.write(f"{r.Chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\t.\t.\n")
    log(f"CADD prescored: {len(hits)} of {len(sites)} indels found, {len(missing)} unscored")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("vep-tables")
    a.add_argument("--raw", required=True, help="split-vep output (fields: " + " ".join(RAW_COLS) + ")")
    a.add_argument("--dbnsfp", required=True, help="slim dbNSFP REVEL/CADD table (bgzipped + tabix)")
    a.add_argument("--popmax", required=True, help="gnomAD popmax > 0.01 variant IDs")
    a.add_argument("--out", required=True, help="output prefix")
    a.add_argument("--all-spliceai", action="store_true",
                   help="score every kept variant with SpliceAI, not just those whose annotation it can change")
    a.set_defaults(func=vep_tables)

    c = sub.add_parser("cadd-prescored")
    c.add_argument("--sites", required=True)
    c.add_argument("--prescored", help="CADD v1.6 GRCh38 prescored indel file (path or URL, tabix-indexed)")
    c.add_argument("--extra", action="append", default=[], help="CADD output file(s) to use first")
    c.add_argument("--out-tsv", required=True)
    c.add_argument("--out-unscored", required=True)
    c.set_defaults(func=cadd_prescored)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
