#!/usr/bin/env python
"""
Stage-by-stage comparison of slim vs original pipeline outputs. Exit code 1 if anything differs.

  compare.py vep      ORIG_TABLE NEW_TABLE        processed VEP tables (transcript rows)
  compare.py spliceai ORIG_VCF   NEW_VCF          max delta score per (variant, gene), on shared variants
  compare.py saige    ORIG_GROUP NEW_GROUP        SAIGE group files
  compare.py cadd     ORIG_TSV   NEW_TSV          CADD PHRED on shared indels
  compare.py self-test                            comparator must flag planted differences

  --no-examples (any position): print counts only, never variant IDs, genes or scores. Use it wherever the
  inputs are derived from individual-level data that must not leave the machine (e.g. UKB on BMRC).
"""

import gzip
import sys
from collections import defaultdict

import pandas as pd

sys.path.insert(0, __file__.rsplit("/validate/", 1)[0] + "/SAIGE_annotations/scripts")
from brava_create_annot import parse_spliceai_info  # noqa: E402

COLS = ["SNP_ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL",
        "BIOTYPE"]


def opener(path):
    return gzip.open(path, "rt", encoding="latin-1") if path.endswith(".gz") else open(path, encoding="latin-1")


def read_table(path):
    df = pd.read_csv(path, sep=" ", dtype=str, keep_default_na=False, encoding="latin-1")
    # --everything also emitted regulatory/motif rows; BRaVa drops them (not protein_coding), the slim run
    # doesn't produce them. Compare transcript rows only.
    df = df[df["TRANSCRIPT"].str.startswith("ENST")]
    df["occ"] = df.groupby(["SNP_ID", "TRANSCRIPT"]).cumcount()
    return df.set_index(["SNP_ID", "TRANSCRIPT", "occ"])


def compare_vep(orig, new, show=5):
    a, b = read_table(orig), read_table(new)
    only_a, only_b = a.index.difference(b.index), b.index.difference(a.index)
    shared = a.index.intersection(b.index)
    print(f"VEP rows: original {len(a)}, slim {len(b)}, shared {len(shared)}, "
          f"only original {len(only_a)}, only slim {len(only_b)}")
    for name, idx in (("only original", only_a), ("only slim", only_b)):
        for k in list(idx[:show]):
            print(f"  {name}: {k}")
    bad = 0
    for c in COLS:
        if c in ("SNP_ID", "TRANSCRIPT"):
            continue
        diff = a.loc[shared, c] != b.loc[shared, c]
        n = int(diff.sum())
        bad += n
        print(f"  {c:12s} mismatches: {n}")
        for k in list(diff[diff].index[:show]):
            print(f"      {k}: original={a.at[k, c]!r} slim={b.at[k, c]!r}")
    return bad + len(only_a) + len(only_b)


def read_spliceai(path):
    out = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
            vid = f"{chrom}:{pos}:{ref}:{alt}"
            fields = [x for x in info.split(";") if x.startswith("SpliceAI=")]
            out[vid] = dict(parse_spliceai_info(fields[0])) if fields else {}
    return out


def compare_spliceai(orig, new, show=5):
    a, b = read_spliceai(orig), read_spliceai(new)
    # chr-prefix agnostic matching
    norm = lambda d: {k if k.startswith("chr") else "chr" + k: v for k, v in d.items()}  # noqa: E731
    a, b = norm(a), norm(b)
    shared = set(a) & set(b)
    print(f"SpliceAI variants: original {len(a)}, slim {len(b)} (slim scores only variants that matter), "
          f"shared {len(shared)}, slim-only {len(set(b) - set(a))}")
    bad = 0
    for v in sorted(shared):
        if a[v] != b[v]:
            bad += 1
            if bad <= show:
                print(f"  {v}: original={a[v]} slim={b[v]}")
    print(f"  variants with any differing gene/max_DS: {bad}")
    return bad + len(set(b) - set(a))


def read_group(path):
    genes = defaultdict(dict)
    with opener(path) as f:
        for line in f:
            gene, kind, *items = line.split()
            genes[gene][kind] = items
    return {g: list(zip(d.get("var", []), d.get("anno", []))) for g, d in genes.items()}


def compare_saige(orig, new, show=5):
    a, b = read_group(orig), read_group(new)
    pairs_a = {(g, v): an for g, l in a.items() for v, an in l}
    pairs_b = {(g, v): an for g, l in b.items() for v, an in l}
    changed = {k for k in set(pairs_a) & set(pairs_b) if pairs_a[k] != pairs_b[k]}
    only_a, only_b = set(pairs_a) - set(pairs_b), set(pairs_b) - set(pairs_a)
    order = sum(1 for g in set(a) & set(b) if a[g] != b[g])
    print(f"SAIGE: genes original {len(a)}, slim {len(b)}; variant-gene pairs original {len(pairs_a)}, "
          f"slim {len(pairs_b)}")
    print(f"  annotation changed: {len(changed)}, only original: {len(only_a)}, only slim: {len(only_b)}, "
          f"genes with any difference (incl. order): {order}")
    trans = defaultdict(int)
    for k in changed:
        trans[(pairs_a[k], pairs_b[k])] += 1
    for (x, y), n in sorted(trans.items(), key=lambda t: -t[1]):
        print(f"    {x} -> {y}: {n}")
    for k in list(sorted(changed))[:show]:
        print(f"    e.g. {k}: {pairs_a[k]} -> {pairs_b[k]}")
    return len(changed) + len(only_a) + len(only_b) + order


def compare_cadd(orig, new, show=5):
    names = ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"]
    a = pd.read_csv(orig, sep="\t", comment="#", header=None, names=names, dtype=str)
    b = pd.read_csv(new, sep="\t", comment="#", header=None, names=names, dtype=str)
    for d in (a, b):
        d["Chrom"] = d["Chrom"].str.replace("^chr", "", regex=True)
    m = a.merge(b, on=["Chrom", "Pos", "Ref", "Alt"], suffixes=("_orig", "_slim"))
    diff = (m["PHRED_orig"].astype(float) - m["PHRED_slim"].astype(float)).abs() > 1e-3
    print(f"CADD indels: original {len(a)}, slim {len(b)}, shared {len(m)}, PHRED differs: {int(diff.sum())}")
    if show and diff.any():
        print(m[diff].head(show).to_string())
    return int(diff.sum())


def self_test():
    """Positive control: each comparator must report planted differences, and none on identical input."""
    import os
    import tempfile
    d = tempfile.mkdtemp()
    table = ("SNP_ID GENE LOF REVEL_SCORE CADD_PHRED CSQ TRANSCRIPT MANE_SELECT CANONICAL BIOTYPE\n"
             "chr1:10:A:G ENSG1 . 0.9 30 missense_variant ENST1 NM_1 YES protein_coding\n"
             "chr1:10:A:G . . . . regulatory_region_variant ENSR1 . . promoter\n")
    group = "ENSG1 var chr1:10:A:G chr1:20:C:T\nENSG1 anno damaging_missense_or_protein_altering synonymous\n"
    vcf = "#h\nchr1\t10\t.\tA\tG\t.\t.\tSpliceAI=G|S---ENSG1.2---ENST1---yes---protein_coding---NM_1|0.10|0.00|0.30|0.00|1|2|3|4\n"
    files = {"t1": table, "t2": table.replace("0.9 30", "0.9 31"), "g1": group,
             "g2": group.replace("synonymous", "non_coding"), "v1": vcf, "v2": vcf.replace("0.30", "0.31")}
    for k, v in files.items():
        with open(os.path.join(d, k), "w") as f:
            f.write(v)
    p = lambda k: os.path.join(d, k)  # noqa: E731
    checks = [("vep identical", compare_vep(p("t1"), p("t1")) == 0),
              ("vep planted", compare_vep(p("t1"), p("t2")) == 1),
              ("saige identical", compare_saige(p("g1"), p("g1")) == 0),
              ("saige planted", compare_saige(p("g1"), p("g2")) > 0),
              ("spliceai identical", compare_spliceai(p("v1"), p("v1")) == 0),
              ("spliceai planted", compare_spliceai(p("v1"), p("v2")) == 1)]
    print()
    for name, ok in checks:
        print(f"self-test {name}: {'PASS' if ok else 'FAIL'}")
    return sum(not ok for _, ok in checks)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--no-examples"]
    show = 5 if len(args) == len(sys.argv) - 1 else 0
    cmd = args[0]
    fn = {"vep": compare_vep, "spliceai": compare_spliceai, "saige": compare_saige, "cadd": compare_cadd}
    n = self_test() if cmd == "self-test" else fn[cmd](*args[1:3], show=show)
    if cmd == "self-test":
        print("SELF-TEST:", "PASSED" if n == 0 else f"{n} FAILED")
    else:
        print("RESULT:", "IDENTICAL" if n == 0 else f"{n} DIFFERENCES")
    sys.exit(1 if n else 0)
