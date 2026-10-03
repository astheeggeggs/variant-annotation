#!/usr/bin/env bash
# Validate the slim pipeline against outputs of the original pipeline for one chromosome.
# Run inside the container (or `pixi run`) on a machine that has the original outputs, e.g. BMRC:
#
#   validate/run_validation.sh -r RESOURCES -o OUT -t 16 [--gpu] \
#     --sites        ORIGINAL_SITES_ONLY_INPUT.vcf.gz          # the VCF that was given to VEP
#     --orig-vep     ..._vep.gnomad_popmax_0.01_processed.txt  # original processed VEP table
#     --orig-spliceai ...sites_only.<chr>.all.vcf             # original SpliceAI output
#     [--orig-cadd   ..._vep_indels.tsv.gz]                    # original CADD indel output, if any
#
# Stages (each prints IDENTICAL or the differences):
#   0. comparator self-test (positive/negative control)
#   1. new vs original brava_create_annot.py on the ORIGINAL intermediates (real-data check of the rewrite)
#   2. slim VEP table vs original (VEP flags, slim cache, dbNSFP table, LOFTEE remote GERP)
#   3. slim SpliceAI vs original on shared variants
#   4. CADD prescored lookup vs original CADD scores
#   5. slim SAIGE file (full SpliceAI, original CADD) vs original-script SAIGE file
#   6. slim SAIGE file in default mode (SpliceAI subset, prescored CADD only) vs the same reference
#   7. (information) the original SAIGE file with SpliceAI cutoff 0.5, and with no SpliceAI at all
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); REPO=$(dirname "$HERE")
THREADS=4; ORIG_CADD=""; GPU=()
while (( $# )); do
  case $1 in
    -r) RES=$2; shift 2 ;; -o) OUT=$2; shift 2 ;; -t) THREADS=$2; shift 2 ;;
    --sites) SITES=$2; shift 2 ;; --orig-vep) OVEP=$2; shift 2 ;;
    --orig-spliceai) OSPL=$2; shift 2 ;; --orig-cadd) ORIG_CADD=$2; shift 2 ;; --gpu) GPU=(--gpu); shift ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
done
mkdir -p "$OUT"; cd "$OUT"
# counts only: the inputs are UKB-derived and nothing identifying may appear in the pasted-back log
cmp_py="python $REPO/validate/compare.py --no-examples"
status=0; run() { echo; echo "### $1"; shift; "$@" || status=1; }

run "0. comparator self-test" $cmp_py self-test

cadd=(); [[ -n $ORIG_CADD ]] && cadd=(--cadd_indels "$ORIG_CADD")
python "$REPO/tests/reference/brava_create_annot_original.py" -v "$OVEP" -s "$OSPL" -w ref.saige.txt "${cadd[@]}" > ref.log 2>&1
python "$REPO/SAIGE_annotations/scripts/brava_create_annot.py" -v "$OVEP" -s "$OSPL" -w new_on_orig.saige.txt "${cadd[@]}" > new_on_orig.log 2>&1
run "1a. rewrite on original intermediates: SAIGE file" cmp ref.saige.txt new_on_orig.saige.txt
run "1b. rewrite on original intermediates: long csv" cmp <(zcat ref.saige.txt.long.csv.gz) <(zcat new_on_orig.saige.txt.long.csv.gz)

extra=(); [[ -n $ORIG_CADD ]] && extra=(--cadd-indels "$ORIG_CADD")
time "$REPO/bin/brava-annotate" -r "$RES" -t "$THREADS" "${GPU[@]}" -o full --keep-work --all-spliceai "${extra[@]}" "$SITES"
name=$(basename "$SITES"); name=${name%.gz}; name=${name%.vcf}; w=full/$name.work

run "2. VEP table (popmax-filtered)" $cmp_py vep "$OVEP" "$w/x.vep.gnomad_popmax_0.01_processed.txt"
run "3. SpliceAI" $cmp_py spliceai "$OSPL" "$w/spliceai.vcf"
if [[ -n $ORIG_CADD ]]; then
  python "$REPO/bin/brava_prep.py" cadd-prescored --sites <(zcat -f "$ORIG_CADD" | awk -v OFS='\t' '!/^#/ {print "chr"$1,$2,$3,$4}') \
    --prescored "$(source "$RES/config.sh"; [[ -s $RES/gnomad.genomes.r3.0.indel.tsv.gz ]] && echo "$RES/gnomad.genomes.r3.0.indel.tsv.gz" \
      || echo https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz)" \
    --out-tsv prescored_all.tsv.gz --out-unscored prescored_missing.vcf
  run "4. CADD prescored vs original CADD" $cmp_py cadd "$ORIG_CADD" prescored_all.tsv.gz
fi
run "5. SAIGE (full SpliceAI, original CADD)" $cmp_py saige ref.saige.txt "full/$name.saige_group.txt"

mkdir -p default/$name.work && cp "$w"/sites.vcf.gz* "$w"/vep.vcf.gz default/$name.work/
time "$REPO/bin/brava-annotate" -r "$RES" -t "$THREADS" "${GPU[@]}" -o default --keep-work "$SITES"
run "6. SAIGE (default: SpliceAI subset, prescored CADD)" $cmp_py saige ref.saige.txt "default/$name.saige_group.txt"
[[ -s default/$name.cadd_unscored_indels.vcf ]] && echo "   (relevant indels without prescored CADD: $(grep -vc '^#' default/$name.cadd_unscored_indels.vcf))"

# Information only (not part of the pass/fail): how much does SpliceAI change the original SAIGE file?
# Same original VEP table and CADD, SpliceAI cutoff 0.5 instead of 0.2, and no SpliceAI at all.
echo; echo "### 7. (information) effect of SpliceAI on the original SAIGE file"
echo "  variant-gene pairs per annotation in the original: $(awk '$2 == "anno" {for (i = 3; i <= NF; i++) n[$i]++}
  END {for (k in n) printf "%s %d; ", k, n[k]}' ref.saige.txt)"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > no_spliceai.vcf
python "$REPO/SAIGE_annotations/scripts/brava_create_annot.py" -v "$OVEP" -s "$OSPL" -w cut05.saige.txt \
  --spliceai_cutoff 0.5 "${cadd[@]}" > cut05.log 2>&1
python "$REPO/SAIGE_annotations/scripts/brava_create_annot.py" -v "$OVEP" -s no_spliceai.vcf -w nosplice.saige.txt \
  "${cadd[@]}" > nosplice.log 2>&1
echo "  -- cutoff 0.2 (original) -> 0.5:"; $cmp_py saige ref.saige.txt cut05.saige.txt | grep -v RESULT || true
echo "  -- cutoff 0.2 (original) -> no SpliceAI:"; $cmp_py saige ref.saige.txt nosplice.saige.txt | grep -v RESULT || true

echo; (( status == 0 )) && echo "ALL STAGES IDENTICAL" || echo "SOME STAGES DIFFER (see above)"
exit $status
