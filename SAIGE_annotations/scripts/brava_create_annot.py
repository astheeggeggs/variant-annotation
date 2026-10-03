#!/usr/bin/env python
# coding: utf-8
"""
Merge VEP with SpliceAI (and optionally CADD indel) information and generate
SAIGE group files according to BRaVa's annotation guidelines.

Drop-in, vectorised replacement for the original row-wise implementation
(kept at tests/reference/brava_create_annot_original.py). Same arguments, same
outputs; tests/test_equivalence.py checks the two agree.
"""

import argparse
import re

import numpy as np
import pandas as pd

PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant",
                 "protein_altering_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

INFRAME_CSQS = ["inframe_deletion", "inframe_insertion"]

CADD_CUTOFF = 28.1
REVEL_CUTOFF = 0.773
SPLICEAI_CUTOFF = 0.20


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge VEP with SpliceAI info and generate annotations according to BRaVa's guidelines.")
    parser.add_argument("--vep", "-v", help="file with VEP annotation (space-delimited)", required=True, type=str)
    parser.add_argument("--spliceai", "-s", help="VCF file with SpliceAI annotations", required=True, type=str)
    parser.add_argument("--out_file", "-w", help="SAIGE output file", required=True, type=str)
    parser.add_argument("--cadd_indels", help="CADD indels file", required=False, type=str)
    parser.add_argument("--spliceai_cutoff", default=SPLICEAI_CUTOFF, type=float,
                        help=f"SpliceAI max delta score at or above which a variant counts as damaging "
                             f"(default {SPLICEAI_CUTOFF}, as in the original)")

    # Columns to read (VEP)
    parser.add_argument("--vep_snp_id_col", default="SNP_ID", help="SNPID (chr:pos:ref:alt) column in VEP table")
    parser.add_argument("--vep_gene_col", default="GENE", help="GENEID column in VEP table")
    parser.add_argument("--vep_lof_col", default="LOF", help="LoF column in VEP table")
    parser.add_argument("--vep_revel_col", default="REVEL_SCORE", help="REVEL column in VEP table")
    parser.add_argument("--vep_cadd_phred_col", default="CADD_PHRED", help="CADD_PHRED column in VEP table")
    parser.add_argument("--vep_consequence_col", default="CSQ", help="Consequence column in VEP table")
    parser.add_argument("--vep_canonical_col", default="CANONICAL", help="Canonical column in VEP table")
    parser.add_argument("--vep_biotype_col", default="BIOTYPE", help="Biotype column in VEP table")
    parser.add_argument("--vep_mane_select_col", default="MANE_SELECT", help="MANE Select column in VEP table")
    return parser.parse_args(argv)


def has_any_csq(csq, terms):
    """True where the '&'-separated consequence string contains any of `terms`."""
    pattern = r"(?:^|&)(?:" + "|".join(map(re.escape, terms)) + r")(?:&|$)"
    return csq.str.contains(pattern, regex=True, na=False)


def get_annotation(df, args):
    """Vectorised version of the original if/elif chain (order matters)."""
    csq = df[args.vep_consequence_col]
    lof = df[args.vep_lof_col]
    missense = has_any_csq(csq, MISSENSE_CSQS)
    inframe = has_any_csq(csq, INFRAME_CSQS)
    other = has_any_csq(csq, OTHER_CSQS)
    synonymous = has_any_csq(csq, SYNONYMOUS_CSQS)
    # NaN comparisons are False, matching the original row-wise behaviour
    damaging_score = (df[args.vep_revel_col] >= REVEL_CUTOFF) | (df[args.vep_cadd_phred_col] >= CADD_CUTOFF)

    conditions = [
        lof == "HC",
        missense & damaging_score,
        df["max_DS"] >= args.spliceai_cutoff,
        lof == "LC",
        missense | inframe,
        other,
        synonymous,
    ]
    choices = [
        "pLoF",
        "damaging_missense_or_protein_altering",
        "damaging_missense_or_protein_altering",
        "damaging_missense_or_protein_altering",
        "other_missense_or_protein_altering",
        "non_coding",
        "synonymous",
    ]
    return pd.Series(np.select(conditions, choices, default=""), index=df.index)


def first_non_missing_revel(value):
    """REVEL comes as '&'-separated per-transcript values; take the first that isn't '.'."""
    if pd.isna(value):
        return np.nan
    for x in str(value).split("&"):
        if x != ".":
            return float(x)
    return np.nan


def parse_spliceai_info(info):
    """
    Return [(ENSG id with version, max DS), ...] for a SpliceAI INFO string, e.g.
    SpliceAI=T|SYMBOL---ENSG0123.4---ENST0123---yes---protein_coding---NM_00123.1|0.00|0.01|0.00|0.03|22|0|-2|8
    Multiple comma-separated gene blocks are handled by the 9-field stride, as in the original.
    """
    fields = info.split("|")
    out = []
    for i in range(len(fields) // 9):
        s = i * 9
        ds_max = max(float(x) for x in fields[s + 2:s + 6])
        out.append((fields[s + 1].split("---")[1], ds_max))
    return out


def read_spliceai(path):
    spliceai_df = pd.read_csv(path, sep="\t", comment="#", header=None,
                              names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
                              na_values=".")
    spliceai_df["ID"] = (spliceai_df.CHROM.astype(str) + ":" + spliceai_df.POS.astype(str) + ":"
                         + spliceai_df.REF + ":" + spliceai_df.ALT)
    spliceai_df = spliceai_df[spliceai_df["INFO"].notna()]

    pairs = spliceai_df["INFO"].map(parse_spliceai_info)
    lengths = pairs.str.len().to_numpy()
    flat = [p for ps in pairs for p in ps]
    exploded = spliceai_df.drop(columns="INFO").loc[spliceai_df.index.repeat(lengths)]
    exploded["GENE"] = [g.split(".")[0] for g, _ in flat]
    exploded["max_DS"] = np.array([d for _, d in flat], dtype=float)

    num_total_rows = exploded.shape[0]
    exploded = exploded[exploded["max_DS"].notna()]
    print(f"DS_Score FILTER: {exploded.shape[0]} out of {num_total_rows} remaining")
    return exploded


def write_saige_file(vep_df, out_file, gene_col, snp_col):
    print(f"Saving the annotation to {out_file}")
    grouped = vep_df.groupby(gene_col, sort=True)
    snps = grouped[snp_col].agg(" ".join)
    annos = grouped["annotation"].agg(" ".join)
    with open(out_file, "w") as fout:
        for gene_id in snps.index:
            fout.write(f"{gene_id} var {snps[gene_id]}\n")
            fout.write(f"{gene_id} anno {annos[gene_id]}\n")


def main(argv=None):
    args = parse_args(argv)

    vep_cols_to_read = [args.vep_snp_id_col, args.vep_gene_col, args.vep_lof_col,
                        args.vep_revel_col, args.vep_cadd_phred_col, args.vep_consequence_col,
                        args.vep_canonical_col, args.vep_biotype_col, args.vep_mane_select_col]

    # cp1252: see https://github.com/BRaVa-genetics/variant-annotation-python/pull/3#issuecomment-1661885882
    vep_df = pd.read_csv(args.vep, sep=" ", usecols=vep_cols_to_read, na_values=".", encoding="cp1252")

    # Keep protein_coding rows on the MANE Select transcript, or the canonical one where MANE is absent
    num_total_rows = vep_df.shape[0]
    mane = vep_df[args.vep_mane_select_col]
    vep_df = vep_df[(vep_df[args.vep_biotype_col] == "protein_coding") &
                    (mane.notna() | (mane.isna() & (vep_df[args.vep_canonical_col] == "YES")))].copy()
    print(f"protein_coding & CANONICAL FILTER: {vep_df.shape[0]} out of {num_total_rows} remaining")

    revel = vep_df[args.vep_revel_col]
    vep_df[args.vep_revel_col] = revel.map(first_non_missing_revel, na_action=None).astype(float)
    print(f"Total REVEL scores: {vep_df[args.vep_revel_col].notna().sum()}")
    print(f"Total pathogenic REVEL scores: {(vep_df[args.vep_revel_col] > REVEL_CUTOFF).sum()}")

    if args.cadd_indels:
        column_names = ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"]
        cadd_indels = pd.read_csv(args.cadd_indels, delimiter="\t", comment="#", names=column_names, header=None)
        cadd_indels["snpid"] = ("chr" + cadd_indels["Chrom"].astype(str) + ":" + cadd_indels["Pos"].astype(str)
                                + ":" + cadd_indels["Ref"] + ":" + cadd_indels["Alt"])
        cadd_scores_before = vep_df[args.vep_cadd_phred_col].notna().sum()
        cadd_indels = cadd_indels.rename(columns={"PHRED": "PHRED_indel"})
        vep_df = pd.merge(vep_df, cadd_indels[["snpid", "PHRED_indel"]],
                          left_on=args.vep_snp_id_col, right_on="snpid", how="left")
        vep_df[args.vep_cadd_phred_col] = vep_df[args.vep_cadd_phred_col].fillna(vep_df["PHRED_indel"])
        cadd_scores_after = vep_df[args.vep_cadd_phred_col].notna().sum()
        print(f"Non-missing CADD PHRED scores: {cadd_scores_before} before, {cadd_scores_after} after indel merge")
        if cadd_scores_before == cadd_scores_after:
            raise Exception("CADD scores unsuccessfully merged")

    spliceai_df = read_spliceai(args.spliceai)
    vep_df = vep_df.merge(spliceai_df, left_on=[args.vep_snp_id_col, args.vep_gene_col],
                          right_on=["ID", "GENE"], how="left")

    # An unsuccessful merge is easy to miss, so fail loudly
    if spliceai_df.shape[0] > 0 and vep_df["max_DS"].notna().sum() == 0:
        raise Exception("""No max_DS SpliceAI scores detected!
            Check variant ID merging between spliceAI and VEP information.
            Perhaps chrCHR:POS:REF:ALT vs CHR:POS:REF:ALT or vice versa""")

    vep_df["annotation"] = get_annotation(vep_df, args)

    num_total_rows = vep_df.shape[0]
    vep_df = vep_df[vep_df["annotation"] != ""]
    print(f"empty annotation FILTER: {vep_df.shape[0]} out of {num_total_rows} remaining")

    print(vep_df["annotation"].value_counts())
    vep_df.to_csv(args.out_file + ".long.csv.gz", index=False, compression="gzip")
    write_saige_file(vep_df, args.out_file, args.vep_gene_col, args.vep_snp_id_col)


if __name__ == "__main__":
    main()
