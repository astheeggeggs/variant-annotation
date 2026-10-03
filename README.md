# BRaVa variant annotation

Turn GRCh38 VCFs into [SAIGE-gene](https://github.com/BRaVa-genetics/universal-saige) group files annotated according to the [BRaVa annotation recommendations](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#), with one command.

```bash
git clone https://github.com/BRaVa-genetics/variant-annotation.git && cd variant-annotation
./run.sh setup -r brava_resources                       # once: ~19 GB (see below)
./run.sh -r brava_resources -t 16 -o out cohort_chr*.vcf.gz
```

`run.sh` uses whatever you have: Apptainer/Singularity, Docker, or [pixi](https://pixi.sh). Nothing else needs installing.

For each input you get `out/<name>.saige_group.txt`, ready for SAIGE-gene, plus:

| file | contents |
| --- | --- |
| `<name>.saige_group.txt.long.csv.gz` | one row per annotated variant–gene, with every score used |
| `<name>.vep_processed.txt.gz` | VEP table for all variants (input to the annotation summaries) |
| `<name>.vep.gnomad_popmax_0.01_processed.txt.gz` | the same table after removing gnomAD popmax > 0.01 variants |
| `<name>.cadd_unscored_indels.vcf` | only written if some indels need CADD scores (see [CADD for indels](#cadd-for-indels)) |

Inputs can be any GRCh38 VCF/BCF: per-chromosome or whole-genome, with or without genotypes, with or without the `chr` prefix, and with multiallelic sites split or not. The pipeline drops genotypes and INFO, splits and left-normalises alleles, and sets IDs to `chr:pos:ref:alt` itself. Interrupted runs resume where they stopped.

## What it does

The steps and versions are the same as the [original pipeline](https://github.com/BRaVa-genetics/variant-annotation/tree/306ab53):

1. **VEP 105 + LOFTEE** (GRCh38 branch) for consequences, MANE Select/canonical transcripts and HC/LC pLoF.
2. **REVEL and CADD (SNVs)** from dbNSFP 4.3a, matched exactly as VEP's dbNSFP plugin matches them.
3. **gnomAD v2.1.1 popmax > 0.01** variants removed.
4. **CADD v1.6 for indels**: only in-frame and protein-altering indels can change an annotation, so only those are looked up.
5. **SpliceAI** with the GENCODE v39/Ensembl 105 annotation. It only scores variants whose annotation SpliceAI can still change; `--all-spliceai` scores everything.
6. `SAIGE_annotations/scripts/brava_create_annot.py` assigns `pLoF`, `damaging_missense_or_protein_altering`, `other_missense_or_protein_altering`, `synonymous` or `non_coding` per gene. It is the same script as before, now vectorised.

### Why it's smaller and faster but gives the same answer

| original | here | why the answer doesn't change |
| --- | --- | --- |
| VEP `--everything` + full 15 GB cache | `--canonical --mane --biotype` + transcript-only cache (~3.5 GB) | The other flags add columns BRaVa never reads (AFs, SIFT, HGVS…) or regulatory/miRNA rows that the protein-coding filter drops. Without a FASTA, VEP was already disabling HGVS. |
| dbNSFP 4.3a (~35 GB, [no longer downloadable](https://dbnsfp.s3.amazonaws.com/dbNSFP4.3a.zip)) via VEP plugin | 2-column table built from the same file | Keeps the first dbNSFP row per (pos, alt, amino-acid change), exactly the row the plugin uses (`pep_match`, ignoring REF), and only for the consequences the plugin annotates. |
| CADD for every indel (220 GB install) | prescored CADD v1.6 lookup for the few indels that can matter | Only missense-class indels that aren't HC pLoF reach a CADD threshold. Any without a prescored score are written out so you can score them with CADD and supply the result. |
| SpliceAI on every variant | SpliceAI on variants whose annotation it can change, in parallel chunks | Variants already called pLoF or damaging missense, and variants with no MANE/canonical protein-coding row, can't change. SpliceAI scores each record independently. |
| popmax filter after VEP | same, applied before SpliceAI/CADD | The filter is by variant ID; nothing upstream depends on other variants. |
| row-wise pandas `apply` | vectorised | `tests/test_equivalence.py` checks byte-identical output against the original script. |

## Resources

`setup` downloads, once:

| resource | size | source |
| --- | --- | --- |
| VEP 105 cache, transcript models only | ~3.5 GB | bundle (or `--vep-cache-from-ensembl`: streams the 15 GB Ensembl tarball and keeps the transcripts) |
| GERP++ bigWig for LOFTEE | 12.6 GB | bundle |
| REVEL/CADD table from dbNSFP 4.3a | ~1 GB | bundle |
| hg38 reference (UCSC) | 1 GB (3 GB unpacked) | UCSC |
| LOFTEE human ancestor + PhyloCSF | 0.9 GB | bundle |

The original needed ~285 GB. On HPC clusters whose compute nodes have no internet, run `setup` on a login node. A shared copy can be reused by everyone via `-r /shared/brava_resources` or `BRAVA_RESOURCES`.

`--bundle URL_OR_DIR` points `setup` at another copy of the bundle, such as a directory on a shared filesystem.

## CADD for indels

SNVs get CADD from dbNSFP. For in-frame and protein-altering indels, `brava-annotate` looks up CADD v1.6 prescored scores remotely (`setup --cadd-prescored local` downloads the 1.2 GB file instead). Any of these indels without a prescored score are written to `<name>.cadd_unscored_indels.vcf`, and the log tells you. Until you supply scores, they are treated as having no CADD score, as in the original pipeline when CADD for indels wasn't run. To match the original exactly, score them with [CADD v1.6](https://github.com/kircherlab/CADD-scripts) (`CADD.sh -g GRCh38 -v v1.6`), or ask the BRaVa team, then rerun with `--cadd-indels scores.tsv.gz` after deleting `<name>.saige_group.txt`.

## GPUs

`./run.sh --gpu ...` uses the GPU image and batched SpliceAI ([this fork](https://github.com/geertvandeweyer/SpliceAI); identical to Illumina's SpliceAI when unbatched). `SPLICEAI_B`/`SPLICEAI_T` set the batch sizes (default 4096/256). On CPU, SpliceAI runs `-t` processes in parallel, each using about 1 GB of RAM.

## Options

```
./run.sh --help
```

## For maintainers

* `tests/test_equivalence.py` (`pixi run -e test test`): new vs original `brava_create_annot.py` on randomised inputs, plus self-tests of the comparators.
* `validate/run_validation.sh`: compares every stage against outputs of the original pipeline for one chromosome (VEP table, SpliceAI, CADD, SAIGE file). Run it on BMRC against the existing UKB outputs before tagging a release.
* `maintainer/build_bundle.sh`: builds the resource bundle from a full install of the original pipeline (VEP cache, dbNSFP 4.3a, LOFTEE data). Upload the output to Zenodo and set `BUNDLE_URL_DEFAULT` in `resources/resources.sh`.
* The container image (`ghcr.io/brava-genetics/variant-annotation:latest[-gpu]`) is built by GitHub Actions from `pixi.lock`, so the container and pixi installs have identical tool versions.

The original step-by-step instructions are in the git history ([`306ab53`](https://github.com/BRaVa-genetics/variant-annotation/tree/306ab53)).
