#!/usr/bin/env bash
# Local stage-2 equivalence test on public data (ClinVar chr21): the slim VEP table vs the original
# pipeline's (vep --everything + full cache + real dbNSFP plugin + LOFTEE, then the upstream
# split-vep/popmax commands). Also exercises maintainer/build_bundle.sh, `brava-annotate setup` and a full
# `brava-annotate` run. Run inside the pixi default env (or the container) with GNU coreutils on PATH.
#
#   validate/local/stage2.sh WORKDIR
#
# WORKDIR must contain (see HANDOFF.md / validate/local/README for how each was made):
#   input.vcf.gz                          test VCF, bare '21' contig names
#   cache_full/homo_sapiens/105_GRCh38/   Ensembl 105 cache: info.txt, chr_synonyms.txt, 21/ (all of it)
#   bundle_src/                           hg38.fa(.fai) [chr21 only], human_ancestor.fa.gz(.fai,.gzi) [21],
#                                         loftee.sql, gerp_conservation_scores.homo_sapiens.GRCh38.bw
#   loftee/                               LOFTEE (loftee_38 at LOFTEE_COMMIT)
#   plugins/dbNSFP.pm                     VEP_plugins release/105
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
W=$(cd "$1" && pwd); cd "$W"
THREADS=${THREADS:-4}
log() { printf '\n[%s] ### %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
TABLE_FMT='%CHROM:%POS:%REF:%ALT %Gene %LoF %REVEL_score %CADD_phred %Consequence %Feature %MANE_SELECT %CANONICAL %BIOTYPE\n'
HEADER='SNP_ID GENE LOF REVEL_SCORE CADD_PHRED CSQ TRANSCRIPT MANE_SELECT CANONICAL BIOTYPE'

# 1. Amino-acid changes for the test variants (any VEP run will do) -> mini dbNSFP with planted edge cases
if [[ ! -s dbnsfp/dbNSFP4.3a_mini.txt.gz.tbi ]]; then
  log "mini dbNSFP"
  mkdir -p dbnsfp
  for c in {1..22} X Y; do printf '%s\tchr%s\n' "$c" "$c"; done > dbnsfp/chr_map.txt
  bcftools annotate --rename-chrs dbnsfp/chr_map.txt input.vcf.gz -Ou \
    | bcftools norm -m-any -f bundle_src/hg38.fa -Oz -o dbnsfp/pre.vcf.gz
  vep -i dbnsfp/pre.vcf.gz -o dbnsfp/pre.vep.vcf --vcf --offline --cache --cache_version 105 --assembly GRCh38 \
    --dir_cache cache_full --no_stats --force_overwrite --quiet --fork "$THREADS"
  bcftools +split-vep dbnsfp/pre.vep.vcf -d -f '%CHROM %POS %REF %ALT %Consequence %Amino_acids\n' \
    | python "$REPO/validate/local/make_mini_dbnsfp.py" > dbnsfp/body.tsv
  (head -1 dbnsfp/body.tsv; tail -n +2 dbnsfp/body.tsv | sort -t$'\t' -k1,1 -k2,2n -s) | bgzip > dbnsfp/dbNSFP4.3a_mini.txt.gz
  tabix -f -s1 -b2 -e2 dbnsfp/dbNSFP4.3a_mini.txt.gz
fi

# 2. Bundle (as a maintainer would on BMRC) and resources (as a user would)
log "bundle"
"$REPO/maintainer/build_bundle.sh" --vep-cache cache_full --dbnsfp dbnsfp/dbNSFP4.3a_mini.txt.gz \
  --loftee-dir bundle_src --out bundle
log "setup"
mkdir -p res && cp -n bundle_src/hg38.fa bundle_src/hg38.fa.fai res/ 2>/dev/null || true  # skip the UCSC download
LOFTEE_PATH="" "$REPO/bin/brava-annotate" setup -r res --bundle bundle
[[ -s res/loftee/LoF.pm ]] || cp -r loftee res/loftee

# 3. Slim
log "slim: brava-annotate"
"$REPO/bin/brava-annotate" -r res -t "$THREADS" -o slim --keep-work input.vcf.gz
SW=slim/input.work

# 4. Original: same sites VCF, upstream VEP command (full cache, --everything, real dbNSFP plugin), upstream
#    post-processing. LOFTEE and dbNSFP.pm share one plugin dir, as in the upstream image.
mkdir -p orig/plugins
cp -r loftee/* orig/plugins/ && cp plugins/dbNSFP.pm orig/plugins/
run_original() {  # $1 = sites VCF, $2 = output prefix
  [[ -s $2.vep.vcf.gz ]] || {
    # upstream's image had BioPerl 1.6.924 (Bio::Perl); the same subs come from lib/perl here
    PERL5LIB="$REPO/lib/perl${PERL5LIB:+:$PERL5LIB}" vep -i "$1" --assembly GRCh38 --vcf --format vcf --cache --dir_cache cache_full -o "$2.vep.vcf" \
      --dir_plugins orig/plugins \
      --plugin "LoF,loftee_path:orig/plugins,human_ancestor_fa:res/human_ancestor.fa.gz,conservation_file:res/loftee.sql,gerp_bigwig:res/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
      --plugin dbNSFP,dbnsfp/dbNSFP4.3a_mini.txt.gz,REVEL_score,CADD_phred \
      --everything --force_overwrite --offline --fork "$THREADS" --no_stats
    bgzip -f "$2.vep.vcf" && tabix -f "$2.vep.vcf.gz"
  }
  bcftools view -i"ID!=@$REPO/resources/gnomad.exomes.r2.1.1.sites.liftover_grch38_popmax_0.01.tsv.bgz" \
    "$2.vep.vcf.gz" -Oz -o "$2.vep.gnomad_popmax_0.01.vcf.gz"
  { echo "$HEADER"; bcftools +split-vep "$2.vep.gnomad_popmax_0.01.vcf.gz" -d -f "$TABLE_FMT"; } \
    > "$2.vep.gnomad_popmax_0.01_processed.txt"
  { echo "$HEADER"; bcftools +split-vep "$2.vep.vcf.gz" -d -f "$TABLE_FMT"; } > "$2.vep.processed.txt"
}
log "original: vep --everything"
run_original "$SW/sites.vcf.gz" orig/chr

status=0
cmp_py() { python "$REPO/validate/compare.py" "$@" || status=1; }
log "compare: popmax-filtered table"
cmp_py vep orig/chr.vep.gnomad_popmax_0.01_processed.txt "$SW/x.vep.gnomad_popmax_0.01_processed.txt"
log "compare: unfiltered table"
cmp_py vep orig/chr.vep.processed.txt "$SW/x.vep_processed.txt"

# 5. Does the original's ANC_ALLELE filter depend on input chr naming? LOFTEE looks up the variant's contig in
#    human_ancestor.fa.gz (bare names). Real ancestral matches are rare, so plant them: set the ancestral base to
#    ALT at every LoF SNV, then run the original on chr21 and on bare 21 input.
log "ANC_ALLELE with chr21 vs 21 input (planted ancestor)"
mkdir -p anc
bcftools +split-vep orig/chr.vep.vcf.gz -d -i 'LoF!="."' -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  | awk 'length($3) == 1 && length($4) == 1' | sort -u > anc/lof_snv.tsv
python - <<'PY'
import gzip
lines = gzip.open("res/human_ancestor.fa.gz", "rt").read().split("\n", 1)
seq = bytearray(lines[1].replace("\n", "").encode())
for l in open("anc/lof_snv.tsv"):
    c, p, r, a = l.split()
    seq[int(p) - 1] = ord(a)
with open("anc/planted.fa", "w") as f:
    f.write(lines[0] + "\n" + "\n".join(seq[i:i + 60].decode() for i in range(0, len(seq), 60)) + "\n")
PY
bgzip -f anc/planted.fa && samtools faidx anc/planted.fa.gz
cut -f1,2 anc/lof_snv.tsv > anc/keep.tsv
bcftools view -T anc/keep.tsv "$SW/sites.vcf.gz" -Oz -o anc/chr.vcf.gz
bcftools annotate --rename-chrs <(printf 'chr21\t21\n') anc/chr.vcf.gz -Oz -o anc/bare.vcf.gz
for n in chr bare; do
  PERL5LIB="$REPO/lib/perl${PERL5LIB:+:$PERL5LIB}" vep -i anc/$n.vcf.gz -o anc/$n.out.vcf --vcf --offline --cache \
    --dir_cache cache_full --assembly GRCh38 --everything --dir_plugins orig/plugins \
    --plugin "LoF,loftee_path:orig/plugins,human_ancestor_fa:anc/planted.fa.gz,conservation_file:res/loftee.sql,gerp_bigwig:res/gerp_conservation_scores.homo_sapiens.GRCh38.bw" \
    --force_overwrite --no_stats --fork "$THREADS" 2> anc/$n.log
  bcftools +split-vep anc/$n.out.vcf -d -f '%ID %Feature %LoF %LoF_filter\n' > anc/$n.lof.txt
  echo "$n input: $(grep -c ANC_ALLELE anc/$n.lof.txt) transcript rows with ANC_ALLELE"
done
if cmp -s anc/chr.lof.txt anc/bare.lof.txt && grep -q ANC_ALLELE anc/chr.lof.txt; then
  echo "ANC_ALLELE: same LoF/LoF_filter for chr21 and 21 input"
else
  echo "ANC_ALLELE: depends on input chr naming (or never fired)"; status=1
fi

echo; (( status == 0 )) && echo "STAGE 2: IDENTICAL" || echo "STAGE 2: DIFFERENCES (see above)"
exit $status
