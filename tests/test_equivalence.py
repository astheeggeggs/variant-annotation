"""
Check the vectorised brava_create_annot.py gives byte-identical outputs to the
original implementation on randomised inputs that exercise the edge cases.

    python tests/test_equivalence.py            # or: pytest tests/
"""

import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ORIGINAL = REPO / "tests/reference/brava_create_annot_original.py"
NEW = REPO / "SAIGE_annotations/scripts/brava_create_annot.py"

TERMS = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
         "frameshift_variant", "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion",
         "missense_variant", "protein_altering_variant", "splice_region_variant", "synonymous_variant",
         "stop_retained_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
         "NMD_transcript_variant", "non_coding_transcript_exon_variant", "upstream_gene_variant",
         "downstream_gene_variant", "coding_sequence_variant", "incomplete_terminal_codon_variant",
         "splice_polypyrimidine_tract_variant"]


def make_inputs(d, seed, n_variants=3000, with_cadd=True):
    rng = random.Random(seed)
    genes = [f"ENSG{n:011d}" for n in range(1, 60)]
    vep_lines = ["SNP_ID GENE LOF REVEL_SCORE CADD_PHRED CSQ TRANSCRIPT MANE_SELECT CANONICAL BIOTYPE"]
    vcf_lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    cadd_lines = ["## CADD GRCh38-v1.6 (c) University of Washington", "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED"]
    pos = 10000
    for _ in range(n_variants):
        pos += rng.randint(1, 500)
        chrom = rng.choice(["chr1", "chr2", "chrX"])
        ref = rng.choice("ACGT")
        kind = rng.random()
        alt = rng.choice([b for b in "ACGT" if b != ref]) if kind < 0.8 else (
            ref + "".join(rng.choice("ACGT") for _ in range(rng.choice([1, 2, 3, 6])))
            if kind < 0.9 else ref)
        if alt == ref:  # deletion
            ref = ref + "".join(rng.choice("ACGT") for _ in range(rng.choice([1, 3, 6])))
        snp = f"{chrom}:{pos}:{ref}:{alt}"
        var_genes = rng.sample(genes, rng.choice([1, 1, 1, 2, 3]))
        for gene in var_genes:
            for t in range(rng.choice([1, 2, 3])):
                csq = "&".join(rng.sample(TERMS, rng.choice([1, 1, 2, 3])))
                lof = rng.choice(["HC", "LC", ".", ".", ".", "."])
                revel = rng.choice([".", ".", f"{rng.random():.3f}",
                                    "&".join(rng.choice([".", f"{rng.random():.3f}"]) for _ in range(rng.randint(1, 4)))])
                cadd = rng.choice([".", f"{rng.uniform(0, 45):.1f}", f"{rng.uniform(25, 35):.3f}"])
                mane = rng.choice([".", ".", f"NM_{rng.randint(1, 99999):06d}.{rng.randint(1, 9)}"])
                canonical = rng.choice(["YES", "."])
                biotype = rng.choice(["protein_coding"] * 4 + ["lncRNA", "nonsense_mediated_decay"])
                vep_lines.append(f"{snp} {gene} {lof} {revel} {cadd} {csq} ENST{rng.randint(1, 10**9):011d} "
                                 f"{mane} {canonical} {biotype}")
        r = rng.random()
        if r < 0.15:
            info = "."
        else:
            blocks = []
            for gene in rng.sample(var_genes + [rng.choice(genes)], rng.choice([1, 1, 2])):
                ds = [rng.choice(["0.00", "0.01", f"{rng.random():.2f}", "0.20", "0.19"]) for _ in range(4)]
                dp = [str(rng.randint(-50, 50)) for _ in range(4)]
                blocks.append(f"{alt}|SYM---{gene}.{rng.randint(1, 20)}---ENST0001---yes---protein_coding---NM_1|"
                              + "|".join(ds + dp))
            info = "SpliceAI=" + ",".join(blocks)
        vcf_lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}")
        if with_cadd and len(ref) != len(alt) and rng.random() < 0.7:
            cadd_lines.append(f"{chrom[3:]}\t{pos}\t{ref}\t{alt}\t{rng.uniform(-1, 8):.6f}\t{rng.uniform(0, 40):.3f}")

    (d / "vep.txt").write_text("\n".join(vep_lines) + "\n")
    (d / "spliceai.vcf").write_text("\n".join(vcf_lines) + "\n")
    if with_cadd:
        (d / "cadd.tsv").write_text("\n".join(cadd_lines) + "\n")


def run(script, d, tag, with_cadd):
    out = d / f"{tag}.saige.txt"
    cmd = [sys.executable, str(script), "-v", str(d / "vep.txt"), "-s", str(d / "spliceai.vcf"), "-w", str(out)]
    if with_cadd:
        cmd += ["--cadd_indels", str(d / "cadd.tsv")]
    subprocess.run(cmd, check=True, capture_output=True)
    with gzip.open(str(out) + ".long.csv.gz", "rt") as f:
        long = f.read()
    return out.read_text(), long


def check(seed, with_cadd):
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        make_inputs(d, seed, with_cadd=with_cadd)
        saige_a, long_a = run(ORIGINAL, d, "orig", with_cadd)
        saige_b, long_b = run(NEW, d, "new", with_cadd)
        assert saige_a == saige_b, f"SAIGE group file differs (seed={seed}, cadd={with_cadd})"
        assert long_a == long_b, f"long csv differs (seed={seed}, cadd={with_cadd})"
        return saige_a.count("\n") // 2, long_a.count("\n") - 1


def test_equivalence():
    for seed in range(5):
        for with_cadd in (True, False):
            check(seed, with_cadd)


if __name__ == "__main__":
    for seed in range(5):
        for with_cadd in (True, False):
            genes, rows = check(seed, with_cadd)
            print(f"seed={seed} cadd={with_cadd}: identical ({genes} genes, {rows} annotated rows)")
