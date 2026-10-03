# Handoff: slim BRaVa variant-annotation (branch `slim`)

Goal: one short command for analysts, minimal downloads/installs, **same SAIGE group files as the original** (upstream `306ab53` + `BRaVa-genetics/vep105_loftee`). Delete this file before merging.

## Decisions taken (with the user)
- CADD for indels: prescored CADD v1.6 lookup by default, plus `--cadd-indels FILE` to supply full CADD scores for the rest (exact).
- Slim resources are built once on BMRC from the existing original install and hosted on Zenodo. Keep dbNSFP **4.3a**: newer releases (v5.x) use GENCODE 46–50, so `pep_match` assignments would change.
- Runtime: one container (Docker/Apptainer) plus a pixi fallback, both from `pixi.lock`.
- GERP: the user chose "remote, full fallback", but **remote turned out not to be viable** (see findings). It is now a local download (12.6 GB).

## Layout
| path | what |
| --- | --- |
| `run.sh` | launcher: picks apptainer/docker/pixi, handles bind mounts, `--gpu` selects the `-gpu` image |
| `bin/brava-annotate` | pipeline driver (setup + run, resumable per step) |
| `bin/brava_prep.py` | rebuilds the original 10-column VEP tables, emulating the dbNSFP plugin from the slim table; picks SpliceAI/CADD targets; CADD prescored lookup |
| `SAIGE_annotations/scripts/brava_create_annot.py` | vectorised drop-in (same CLI/outputs) |
| `resources/resources.sh` | resource URLs, `setup`, checks. **`BUNDLE_URL_DEFAULT` is a TBD placeholder** |
| `resources/*.gz/.bgz` | gnomAD popmax list (moved from vep105_loftee) and GENCODE v39 SpliceAI annotation (moved from `data/SpliceAI`) |
| `maintainer/build_bundle.sh` | builds the bundle (transcript-only VEP cache tar, dbNSFP→REVEL/CADD table, LOFTEE files, MD5SUMS) |
| `validate/compare.py`, `validate/run_validation.sh` | stage-by-stage comparison against original outputs on BMRC (with comparator self-test) |
| `validate/local/` | local stage-2 test on public data: `stage2.sh WORKDIR` (driver), mini dbNSFP generator, synthetic GERP bigWig, bigWig query test |
| `lib/perl/Bio/Perl.pm` | `reverse_complement` from BioPerl 1.6.924 for LOFTEE (conda BioPerl 1.7.8 lacks Bio::Perl) |
| `bin/spliceai_cpu.py` | CPU SpliceAI: the fork's record-by-record loop (the fork's own `-B 1` path crashes) |
| `tests/` | original script kept as oracle + randomised equivalence test |
| `Dockerfile`, `pixi.toml`, `pixi.lock`, `.github/workflows/docker.yml` | env/image; CI runs tests and publishes `ghcr.io/<repo>:latest[-gpu]` |

## Verified
- `brava_create_annot.py` rewrite: byte-identical SAIGE file and long.csv vs the original on 10 randomised datasets. Negative controls (cutoff 0.20→0.19; HC/damaging precedence swap) are detected. `pixi run -e test test`.
- `validate/compare.py self-test`: no false positives, and every planted difference is flagged.
- `pixi.lock` solves (linux-64, osx-64; gpu env linux-64 only).
- Docker image builds (3.6 GB uncompressed). Inside it: VEP 105.0, bcftools/htslib/samtools 1.21, split-vep, LOFTEE perl modules, SpliceAI import with TF/Keras 2.15.

- **Local stage-2 equivalence (2026-10-03)**, `validate/local/stage2.sh` on 12,107 ClinVar chr21 variants (bare `21` names, 28 made multiallelic) with a mini dbNSFP carrying planted edge cases (613 later duplicates, 543 `.`-first rows, 700 decoy aa changes, 379 wrong REF). The slim VEP table is **identical** to the original (`vep --everything`, full cache, real release/105 dbNSFP plugin, LOFTEE, upstream popmax/split-vep): 95,929 rows unfiltered and 93,260 after popmax, 0 mismatches in every column. The same run built the bundle with `build_bundle.sh` and ran `brava-annotate setup --bundle DIR`.
- **ANC_ALLELE does not depend on input chr naming.** With the ancestral base planted as ALT at all 452 LoF SNVs, the original gives identical LoF/LoF_filter for `chr21` and `21` input (1,712 ANC_ALLELE rows each). VEP passes the cache's contig name to LOFTEE. No real LoF SNV in the sample had an ancestral ALT.
- Chunked parallel SpliceAI output is byte-identical to a single process (200 variants, 4 chunks).
- **Local end to end**: `brava-annotate` (default mode, CPU SpliceAI on 8,775 variants in 6 chunks) gives a SAIGE group file and long csv **byte-identical** to the original script run on the original pipeline's VEP table with the same SpliceAI/CADD inputs (210 genes, 13,456 variant-gene pairs).

## NOT yet verified (in priority order)
1. **Docker image**: rebuild and re-check that LOFTEE runs in it (before the Bio::Perl fix it silently did not). `brava-annotate` itself has now run end to end locally (pixi, macOS, bash 5 + GNU coreutils), CPU SpliceAI only; the `--gpu` path is untested.
2. ~~Local stage-2 equivalence test~~ done, see Verified. Setup notes kept for reruns: Needs ~25 GB free plus the chr21 VEP cache:
   - Stream the Ensembl 105 cache and extract only `homo_sapiens/105_GRCh38/{21/*,info.txt,chr_synonyms.txt}` (chr21 arrives late in the tar; ~2 h at 2 MB/s).
   - Test resources: UCSC `chr21.fa` as `hg38.fa`; chr21 of `human_ancestor` (`samtools faidx <https URL> 21`, renamed `21`); `loftee.sql`; a synthetic GERP bigWig for contig `21` (pyBigWig, bioconda). The GERP/ancestor files only need to be identical between the two runs.
   - Input: ClinVar GRCh38 chr21 (`bcftools view -r 21 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`), subsampled (~12k).
   - Run a slim VEP once to get amino-acid changes → `validate/local/make_mini_dbnsfp.py` → `dbNSFP4.3a_mini.txt.gz` (the filename must match `/4\./`) → table via the awk in `maintainer/build_bundle.sh`.
   - "Original": `vep --everything` + full chr21 cache + real dbNSFP plugin (`VEP_plugins` release/105 `dbNSFP.pm`) on the mini file + LOFTEE → upstream split-vep/popmax commands. "Slim": `brava-annotate --keep-work`. Compare with `validate/compare.py vep`.
   - Also check that chunked parallel SpliceAI output is byte-identical to a single process on a few hundred variants.
3. **Real-data validation on BMRC**: `validate/run_validation.sh` against the existing UKB chr21 outputs of the original pipeline (VEP processed table, SpliceAI VCF, CADD indel tsv).
4. Build the bundle on BMRC (`maintainer/build_bundle.sh`), upload to Zenodo, and set `BUNDLE_URL_DEFAULT`. Check that REVEL/CADD/dbNSFP non-commercial terms allow redistributing the derived table.

## Findings to carry forward
- **SpliceAI writes `.` scores** for variants it doesn't score (REF and ALT both > 1 bp). The original script crashes on them (it meant to drop them); the rewrite treats them as missing. It also uses `np.max` like the original, because SpliceAI writes `-0.00` and the sign of a zero maximum shows in the long csv.
- **SpliceAI's effect, local ClinVar chr21** (REVEL/CADD synthetic here, so only SpliceAI/LOFTEE/ClinVar are real): SpliceAI promotes 380 of 13,456 variant-gene pairs (98 of 210 genes) to damaging at 0.2, and 146 at 0.5. Among ClinVar-classified variants, P/LP of P+B: non-coding fallback DS >= 0.5: 64% (14/22) vs 0.2-0.5: 9% (3/32) vs < 0.2: 1%; synonymous fallback >= 0.5: 20% (2/10), 0.2-0.5: 6% (2/35); HC pLoF benchmark 98%. Suggests >= 0.5 carries signal and most 0.2-0.5 promotions look benign. Small n and ClinVar ascertainment; stage 7 of `run_validation.sh` gives the UKB counts.
- **LOFTEE needs `Bio::Perl`, which the conda BioPerl 1.7.8 packages don't ship.** VEP then only *warns* ("Failed to compile plugin LoF") and writes output with no LoF field. `run_vep` now prepends `lib/perl` to PERL5LIB and dies if the CSQ header lacks `LoF`.
- **SpliceAI fork `-B 1` path crashes** (`args_output_data` typo in `run_spliceai`, present since c7fd579; the pinned rev is the fork's HEAD). CPU runs go through `bin/spliceai_cpu.py`, which uses the same `spliceai.utils` scoring. The fork's `get_delta_scores` is a refactor of Illumina 1.3.1 with the same arithmetic.
- The VEP dbNSFP plugin requires an `Ensembl_transcriptid` column in the header (fixed in `make_mini_dbnsfp.py`).
- `maintainer/build_bundle.sh` was committed without the executable bit (fixed).
- Under Rosetta, x86 TensorFlow SpliceAI scores about 8 variants/min per process, after a model load of about 2.5 min.
- **dbNSFP 4.3a download URL in upstream `download_data.sh` is 404.** New cohorts can't reproduce the original, which is what motivates the hosted table.
- **Remote GERP is not viable.** `Bio::DB::BigWig` only accepts `http:`/`ftp:`, and the bioconda Kent lib has no OpenSSL, so it can't read Broad or Zenodo (https). Plain http works (tested against UCSC), but LOFTEE reopens the file for every exon interval at ~0.5 s per open, so a chromosome would take days. Follow-up idea for an *exact* small GERP: a copy with the same byte layout where data blocks outside protein-coding CDS are replaced by valid empty blocks. END_TRUNC only uses GERP when the variant is ≤50 bp from the last exon. Needs care: LOFTEE queries GERP for every stop-gained/frameshift transcript, and a corrupt block aborts. The file has zoom levels from 40 bp, so a naively re-written subset is *not* exact.
- VEP 105 with `--everything --offline` and no FASTA disables HGVS ("INFO: Disabling --hgvs"). The original therefore ran LOFTEE without a FASTA, so the slim run must **not** pass `--fasta`. `build_bundle.sh` also excludes any FASTA from the cache tar, because VEP auto-loads one.
- dbNSFP plugin (VEP 105): annotates only missense/stop_gained/stop_lost/start_lost TVAs and SNVs only. It matches `pos`, `alt` and `aaref/aaalt` (X→*) against `Amino_acids`, ignores REF, takes the **first** matching row, and drops `.` values. VEP's VCF writer turns `;`/`,`/`|` into `&`. The python takes the first non-`.` REVEL.
- VEP cache composition (chrs 1, 10–18 measured): transcript 1.15 GB, regulatory 1.04 GB, variation 3.71 GB. Transcript-only for the whole genome is ~3–3.5 GB of 15.4 GB. Further trim possible: strip SIFT/PolyPhen matrices from the transcript objects (needs validation).
- Upstream image `ghcr.io/brava-genetics/vep105_loftee:main` **is** publicly pullable. The 403 seen locally came from the Docker Desktop daemon on this machine (it hit every ghcr image; anonymous curl got tokens fine). For the same reason the Dockerfile uses `ubuntu:24.04` + the pixi install script rather than a ghcr base.
- `perl-bioperl` 1.7.8 (needed by bioconda ensembl-vep 105) drags in openjdk/gcc/mysql/blast. Pinning 1.7.2 conflicts with perl 5.32. The Dockerfile prunes compilers, Java, docs, headers and static libs (image 4.85 → 3.6 GB). Re-check that the pruned image still runs the full pipeline.
- `spliceai` needs `setuptools<81` (pkg_resources).
- Speeds measured from Oxford: Ensembl ~2 MB/s, Zenodo ~0.4 MB/s, UCSC ~60 KB/s, `personal.broadinstitute.org` ~80 KB/s with frequent resets. Mirror the LOFTEE files in the bundle.
- `pixi run '…'` re-parses quotes (deno shell); environments with `no-default-feature` don't see top-level tasks (hence `[feature.test.tasks]`).

## Known gaps / TODO in code
- `run.sh` uses `mapfile` (bash ≥4). macOS system bash 3.2 will fail; fine on Linux/HPC.
- In default mode, `long.csv.gz` has NaN `max_DS` for variants SpliceAI was skipped on (already pLoF/damaging), and NaN CADD for irrelevant indels. The SAIGE file is unaffected; `--all-spliceai` restores the fuller column.
- `validate/run_validation.sh` stage 6 is *expected* to differ only by in-frame indels lacking prescored CADD. It reports the transitions.
- README validation claims must be updated once items 1–3 pass.
