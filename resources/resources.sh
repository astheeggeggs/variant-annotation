# Resource locations, one-off setup and checks for brava-annotate (sourced, not executed).
#
# Everything except the reference genome comes from one "bundle" (a Zenodo record, or a local
# directory such as a shared HPC path) built by maintainer/build_bundle.sh. The bundle replaces
# ~285 GB of upstream downloads with ~19 GB (12.6 GB of which is the GERP bigWig):
#   vep105_GRCh38_transcripts.tar          VEP 105 cache, transcript models only (no variation/regulation)
#   dbnsfp4.3a_revel_cadd.tsv.gz(.tbi)     REVEL_score + CADD_phred from dbNSFP 4.3a, first row per plugin key
#   human_ancestor.fa.gz(.fai,.gzi), loftee.sql.gz, gerp_conservation_scores.homo_sapiens.GRCh38.bw   LOFTEE
#   MD5SUMS

# TODO(maintainers): replace with the Zenodo record once the bundle is uploaded
BUNDLE_URL_DEFAULT=${BRAVA_BUNDLE:-"https://zenodo.org/records/TBD/files"}

HG38_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
ENSEMBL_CACHE_URL="https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz"
LOFTEE_GIT="https://github.com/populationgenomics/loftee_38.git"
LOFTEE_COMMIT="c8fdde00e515148450416128d43fcf01f1ee6bb8"
GERP_FILE="gerp_conservation_scores.homo_sapiens.GRCh38.bw"
CADD_PRESCORED_FILE="gnomad.genomes.r3.0.indel.tsv.gz"
CADD_PRESCORED_URL="https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/$CADD_PRESCORED_FILE"

VEP_CACHE_TAR="vep105_GRCh38_transcripts.tar"
DBNSFP_TABLE="dbnsfp4.3a_revel_cadd.tsv.gz"
POPMAX_LIST="gnomad.exomes.r2.1.1.sites.liftover_grch38_popmax_0.01.tsv.bgz"   # shipped in the repo
GENCODE_ANNOT="gencode.v39.ensembl.v105.annotation.txt.gz"                       # shipped in the repo

# Fetch $1 (relative to the bundle, or an absolute URL) to $2
fetch() {
  local src=$1 dest=$2
  [[ $src == http* || $src == /* ]] || src=$BUNDLE/$src
  if [[ $src == http* ]]; then
    [[ $src == *TBD* ]] && die "the resource bundle isn't hosted yet; pass --bundle URL_OR_DIR (see README)"
    log "  downloading $src"
    curl -fL --retry 10 --retry-all-errors -C - -o "$dest.part" "$src" && mv "$dest.part" "$dest"
  else
    [[ -e $src ]] || die "missing $src"
    log "  copying $src"
    cp "$src" "$dest"
  fi
}

md5_check() {
  local f=$1 sums=$RESOURCES/MD5SUMS expected
  [[ -s $sums ]] || return 0
  expected=$(awk -v f="$(basename "$f")" '$2==f {print $1}' "$sums")
  [[ -n $expected ]] || return 0
  [[ $(md5sum "$f" | cut -d' ' -f1) == "$expected" ]] || die "checksum mismatch for $f (delete it and rerun setup)"
}

setup_resources() {
  local cadd_mode=remote vep_from_ensembl=0
  BUNDLE=$BUNDLE_URL_DEFAULT
  while (( $# )); do
    case $1 in
      -r|--resources) RESOURCES=$2; shift 2 ;;
      --bundle) BUNDLE=$2; shift 2 ;;
      --cadd-prescored) cadd_mode=$2; shift 2 ;;
      --vep-cache-from-ensembl) vep_from_ensembl=1; shift ;;
      -h|--help) usage; return 0 ;;
      *) die "unknown setup option $1" ;;
    esac
  done
  [[ $BUNDLE == http* ]] || BUNDLE=$(cd "$BUNDLE" && pwd)
  mkdir -p "$RESOURCES"; RESOURCES=$(cd "$RESOURCES" && pwd)
  log "Setting up resources in $RESOURCES (bundle: $BUNDLE)"
  cd "$RESOURCES"

  [[ -s MD5SUMS ]] || fetch MD5SUMS MD5SUMS || true

  # Reference genome (bcftools norm + SpliceAI). SpliceAI's pyfaidx needs it uncompressed.
  if [[ ! -s hg38.fa.fai ]]; then
    fetch "$HG38_URL" hg38.fa.gz
    gunzip -f hg38.fa.gz && samtools faidx hg38.fa
  fi

  # VEP cache: transcript models only
  if [[ ! -s vep/homo_sapiens/105_GRCh38/info.txt ]]; then
    mkdir -p vep
    if (( vep_from_ensembl )); then
      log "  streaming the full Ensembl cache (15 GB) and keeping transcript models only"
      curl -fL --retry 10 "$ENSEMBL_CACHE_URL" \
        | tar -xzf - -C vep --exclude='*/all_vars.gz*' --exclude='*_reg.gz' --exclude='*_var.gz'
    else
      fetch "$VEP_CACHE_TAR" "$VEP_CACHE_TAR" && md5_check "$VEP_CACHE_TAR"
      tar -xf "$VEP_CACHE_TAR" -C vep && rm -f "$VEP_CACHE_TAR"
    fi
  fi

  local f
  for f in "$DBNSFP_TABLE" "$DBNSFP_TABLE.tbi" human_ancestor.fa.gz human_ancestor.fa.gz.fai human_ancestor.fa.gz.gzi; do
    [[ -s $f ]] || { fetch "$f" "$f"; md5_check "$f"; }
  done
  if [[ ! -s loftee.sql ]]; then fetch loftee.sql.gz loftee.sql.gz && md5_check loftee.sql.gz && gunzip -f loftee.sql.gz; fi

  # GERP bigWig for LOFTEE's END_TRUNC filter (12.6 GB). Must be local: LOFTEE reopens the file for every
  # exon interval, and the Kent library in perl-bio-bigfile has no https support, so remote reads are
  # impractically slow (~0.5 s per open over plain http) or impossible (https).
  [[ -s $GERP_FILE ]] || { fetch "$GERP_FILE" "$GERP_FILE"; md5_check "$GERP_FILE"; }
  if [[ $cadd_mode == local ]]; then
    for f in "$CADD_PRESCORED_FILE" "$CADD_PRESCORED_FILE.tbi"; do [[ -s $f ]] || fetch "$CADD_PRESCORED_URL${f#$CADD_PRESCORED_FILE}" "$f"; done
  fi

  # LOFTEE plugin code (already in the container image at /opt/loftee)
  if [[ ! -d /opt/loftee && ! -s loftee/LoF.pm ]]; then
    git clone -q "$LOFTEE_GIT" loftee && git -C loftee checkout -q "$LOFTEE_COMMIT"
  fi

  cat > config.sh <<EOF
CADD_MODE=$cadd_mode
BUNDLE=$BUNDLE
EOF
  touch .complete
  log "Resources ready: $(du -sh "$RESOURCES" | cut -f1) in $RESOURCES"
}

check_resources() {
  [[ -e $RESOURCES/.complete ]] || die "resources not set up in $RESOURCES; run: brava-annotate setup -r $RESOURCES"
  RESOURCES=$(cd "$RESOURCES" && pwd)
  # shellcheck source=/dev/null
  source "$RESOURCES/config.sh"
  command -v vep >/dev/null || die "vep not found; use the container or 'pixi run' (see README)"
}

gerp_path() { echo "$RESOURCES/$GERP_FILE"; }

cadd_prescored_path() {
  if [[ -s $RESOURCES/$CADD_PRESCORED_FILE ]]; then echo "$RESOURCES/$CADD_PRESCORED_FILE"
  else echo "$CADD_PRESCORED_URL"; fi
}
