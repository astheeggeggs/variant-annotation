#!/usr/bin/env bash
# Build the slim resource bundle from an existing full install of the original pipeline
# (e.g. BMRC's vep105_loftee/vep_data). Run once by maintainers; upload the output to Zenodo.
#
#   maintainer/build_bundle.sh \
#     --vep-cache  /path/to/vep_data            # contains homo_sapiens/105_GRCh38/
#     --dbnsfp     /path/to/dbNSFP4.3a.txt.gz   # as made by vep105_loftee/download_data.sh
#     --loftee-dir /path/to/vep_data            # human_ancestor.fa.gz(.fai,.gzi), loftee.sql(.gz), gerp bigWig
#     --out        bundle/
set -euo pipefail

THREADS=4
while (( $# )); do
  case $1 in
    --vep-cache) VEP_CACHE=$2; shift 2 ;;
    --dbnsfp) DBNSFP=$2; shift 2 ;;
    --loftee-dir) LOFTEE_DIR=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac
done
: "${VEP_CACHE:?--vep-cache required}" "${DBNSFP:?--dbnsfp required}" "${LOFTEE_DIR:?--loftee-dir required}" "${OUT:?--out required}"
mkdir -p "$OUT"
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# 1. VEP cache: transcript models only. Variation (all_vars) and regulatory (_reg) data are only read
#    with --check_existing/--af*/--regulatory, none of which affect the columns BRaVa uses. Any FASTA
#    in the cache dir is excluded too: VEP would auto-load it, and the original run had none.
if [[ ! -s $OUT/vep105_GRCh38_transcripts.tar ]]; then
  log "VEP cache -> transcript-only tar"
  tar -cf "$OUT/vep105_GRCh38_transcripts.tar.part" -C "$VEP_CACHE" \
    --exclude='all_vars.gz*' --exclude='*_reg.gz' --exclude='*_var.gz' \
    --exclude='*.fa' --exclude='*.fa.gz*' --exclude='*.fa.fai' homo_sapiens/105_GRCh38
  mv "$OUT/vep105_GRCh38_transcripts.tar.part" "$OUT/vep105_GRCh38_transcripts.tar"
fi

# 2. dbNSFP -> REVEL/CADD table that reproduces the VEP 105 dbNSFP plugin (pep_match=1, default
#    consequences). The plugin uses the FIRST row matching (pos, alt, aaref/aaalt with X->*), ignoring
#    REF, and drops '.' values; VEP then writes ';' and '|' as '&'. We keep exactly that row per key and
#    drop keys where both values are '.', which the plugin would also emit as missing.
if [[ ! -s $OUT/dbnsfp4.3a_revel_cadd.tsv.gz.tbi ]]; then
  log "dbNSFP -> REVEL/CADD table (streams the whole file; ~1-2 h)"
  zcat -f "$DBNSFP" | awk -F'\t' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      split("#chr pos(1-based) ref alt aaref aaalt REVEL_score CADD_phred", need, " ")
      for (j in need) if (!(need[j] in col)) { print "missing column " need[j] > "/dev/stderr"; exit 1 }
      c = col["#chr"]; p = col["pos(1-based)"]; r = col["ref"]; a = col["alt"]
      ar = col["aaref"]; aa = col["aaalt"]; rv = col["REVEL_score"]; cd = col["CADD_phred"]
      print "#chr", "pos", "ref", "alt", "aa", "REVEL_score", "CADD_phred"
      next
    }
    {
      if ($c != last_chr || $p != last_pos) { delete seen; last_chr = $c; last_pos = $p }
      pep = $ar "/" $aa; gsub(/X/, "*", pep)
      k = $a SUBSEP pep
      if (k in seen) next
      seen[k] = 1
      rev = $rv; cad = $cd
      if (rev == "." && cad == ".") next
      gsub(/[;|,]/, "\\&", rev); gsub(/[;|,]/, "\\&", cad)
      gsub(/[ \t]+/, "_", rev); gsub(/[ \t]+/, "_", cad)
      print $c, $p, $r, $a, pep, rev, cad
    }' | bgzip -@ "$THREADS" > "$OUT/dbnsfp4.3a_revel_cadd.tsv.gz.part"
  mv "$OUT/dbnsfp4.3a_revel_cadd.tsv.gz.part" "$OUT/dbnsfp4.3a_revel_cadd.tsv.gz"
  tabix -f -s1 -b2 -e2 "$OUT/dbnsfp4.3a_revel_cadd.tsv.gz"
fi

# 3. LOFTEE data, copied as-is (mirrored so we don't depend on a personal web page)
for f in human_ancestor.fa.gz human_ancestor.fa.gz.fai human_ancestor.fa.gz.gzi gerp_conservation_scores.homo_sapiens.GRCh38.bw; do
  [[ -s $OUT/$f ]] || { log "copy $f"; cp "$LOFTEE_DIR/$f" "$OUT/$f"; }
done
if [[ ! -s $OUT/loftee.sql.gz ]]; then
  if [[ -s $LOFTEE_DIR/loftee.sql.gz ]]; then cp "$LOFTEE_DIR/loftee.sql.gz" "$OUT/"; else gzip -c "$LOFTEE_DIR/loftee.sql" > "$OUT/loftee.sql.gz"; fi
fi

log "checksums"
( cd "$OUT" && md5sum vep105_GRCh38_transcripts.tar dbnsfp4.3a_revel_cadd.tsv.gz dbnsfp4.3a_revel_cadd.tsv.gz.tbi \
    human_ancestor.fa.gz human_ancestor.fa.gz.fai human_ancestor.fa.gz.gzi loftee.sql.gz \
    gerp_conservation_scores.homo_sapiens.GRCh38.bw > MD5SUMS )
du -ch "$OUT"/* | sort -h
log "Done. Upload everything in $OUT to one Zenodo record, then set BUNDLE_URL_DEFAULT in resources/resources.sh"
