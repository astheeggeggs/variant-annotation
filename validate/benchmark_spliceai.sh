#!/usr/bin/env bash
# SpliceAI speed on this machine, CPU and (optionally) GPU, and whether GPU scores match CPU scores.
# Uses random SNVs in GENCODE v39 exons +-50 bp generated from the reference (no downloads, no cohort data),
# so the output is safe to paste anywhere. Run in the pixi env (`-e gpu` for --gpu) or the container:
#
#   pixi run validate/benchmark_spliceai.sh -R /path/hg38.fa -t 16                  # CPU node
#   pixi run -e gpu validate/benchmark_spliceai.sh -R /path/hg38.fa -t 4 --gpu      # GPU node
#
# GPU batch sizes as in brava-annotate: SPLICEAI_B (default 4096) and SPLICEAI_T (default 256, sized for a
# 40 GB card; try 64 on a 10-12 GB GPU slice if it runs out of memory).
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SPLICEAI_B=${SPLICEAI_B:-4096}; SPLICEAI_T=${SPLICEAI_T:-256}
THREADS=4; GPU=0; N_GPU=20000; N_CPU=200; OUT=spliceai_benchmark
while (( $# )); do
  case $1 in
    -R) REF=$2; shift 2 ;; -t) THREADS=$2; shift 2 ;; --gpu) GPU=1; shift ;;
    -n) N_GPU=$2; shift 2 ;; -o) OUT=$2; shift 2 ;;
    *) echo "usage: $0 -R hg38.fa [-t THREADS] [--gpu] [-n GPU_VARIANTS] [-o OUTDIR]"; exit 1 ;;
  esac
done
ANN=$REPO/resources/gencode.v39.ensembl.v105.annotation.txt.gz
mkdir -p "$OUT"; cd "$OUT"

# Random SNVs in exons +-50 bp (the region a precomputed table would cover), sorted, fixed seed
python - "$REF" "$ANN" "$(( N_GPU > N_CPU * THREADS ? N_GPU : N_CPU * THREADS ))" <<'PY'
import bisect, gzip, random, sys
import pysam
ref, ann, n = pysam.FastaFile(sys.argv[1]), sys.argv[2], int(sys.argv[3])
iv = []
with gzip.open(ann, "rt") as f:
    next(f)
    for l in f:
        x = l.split("\t")
        if x[1] in ref.references:
            iv += [(x[1], int(s) - 50, int(e) + 50) for s, e in zip(x[5].strip(",").split(","), x[6].strip(",").split(","))]
cum = [0]
for c, s, e in iv:
    cum.append(cum[-1] + e - s)
rng, out = random.Random(1), set()
while len(out) < n:
    c, s, e = iv[bisect.bisect_right(cum, rng.randrange(cum[-1])) - 1]
    pos = rng.randrange(s, e) + 1
    b = ref.fetch(c, pos - 1, pos).upper()
    if b in "ACGT":
        out.add((c, pos, b, rng.choice([a for a in "ACGT" if a != b])))
key = lambda v: (ref.references.index(v[0]), v[1], v[3])
# random order, so any first-k subset is a random sample (SpliceAI scores records independently)
vs = sorted(out, key=key)
rng.shuffle(vs)
with open("all.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in sorted({v[0] for v in vs}, key=ref.references.index):
        f.write(f"##contig=<ID={c},length={ref.get_reference_length(c)}>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for c, p, r, a in vs:
        f.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n")
PY

score_cpu() { CUDA_VISIBLE_DEVICES="" TF_NUM_INTRAOP_THREADS=1 TF_NUM_INTEROP_THREADS=1 OMP_NUM_THREADS=1 \
  python "$REPO/bin/spliceai_cpu.py" -I "$1" -O "$2" -R "$REF" -A "$ANN" 2> "$2.log"; }
subset() { { grep '^#' all.vcf; grep -v '^#' all.vcf | sed -n "$1"; } > "$2"; }
n_pred() { grep -o 'SpliceAI=[^[:space:]]*' "$1" | tr ',' '\n' | grep -c '|' || true; }
report() { awk -v m="$1" -v n="$2" -v p="$3" -v t="$4" 'BEGIN {
  printf "%-24s %6d variants %7d predictions %6ds -> %9.0f predictions/hour\n", m, n, p, t, p * 3600 / (t ? t : 1) }'; }
echo "host: $(hostname -s), CPUs: $(nproc), GPU: $( (( GPU )) && (nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) || echo none)"

# CPU, one process: 10 variants (~ startup cost) and N_CPU variants
subset "1,10p" c10.vcf; subset "1,${N_CPU}p" c1.vcf
t=$SECONDS; score_cpu c10.vcf c10.out.vcf; t10=$(( SECONDS - t ))
t=$SECONDS; score_cpu c1.vcf c1.out.vcf; t1=$(( SECONDS - t ))
echo "CPU, 10 variants (mostly start-up: loading the 5 models): ${t10}s"
report "CPU, 1 process" "$N_CPU" "$(n_pred c1.out.vcf)" "$t1"

# CPU, THREADS processes in parallel (what brava-annotate does), N_CPU variants each
if (( THREADS > 1 )); then
  t=$SECONDS
  for i in $(seq 0 $(( THREADS - 1 ))); do
    subset "$(( i * N_CPU + 1 )),$(( (i + 1) * N_CPU ))p" "p$i.vcf"; score_cpu "p$i.vcf" "p$i.out.vcf" &
  done; wait
  tp=$(( SECONDS - t )); np=0; for i in $(seq 0 $(( THREADS - 1 ))); do np=$(( np + $(n_pred "p$i.out.vcf") )); done
  report "CPU, $THREADS processes" "$(( N_CPU * THREADS ))" "$np" "$tp"
fi

if (( GPU )); then
  subset "1,${N_GPU}p" g.vcf
  t=$SECONDS
  spliceai -I g.vcf -O g.out.vcf -R "$REF" -A "$ANN" -B "$SPLICEAI_B" -T "$SPLICEAI_T" -t "$PWD" > g.log 2>&1 \
    || { tail -5 g.log; echo "GPU SpliceAI failed (log: $PWD/g.log)"; exit 1; }
  tg=$(( SECONDS - t ))
  report "GPU (-B $SPLICEAI_B -T $SPLICEAI_T)" "$N_GPU" "$(n_pred g.out.vcf)" "$tg"
  # Do GPU and CPU agree on the variants both scored?
  python - c1.out.vcf g.out.vcf <<'PY'
import sys
def read(p):
    d = {}
    for l in open(p):
        if not l.startswith("#"):
            x = l.rstrip("\n").split("\t")
            s = [f for f in x[7].split(";") if f.startswith("SpliceAI=")]
            d[tuple(x[:2] + x[3:5])] = s[0] if s else ""
    return d
a, b = read(sys.argv[1]), read(sys.argv[2])
shared = set(a) & set(b)
diff = [k for k in shared if a[k] != b[k]]
def ds(s):
    return [float(v) for e in s[9:].split(",") for v in e.split("|")[2:6] if v != "."] if s else []
dmax = max((max(abs(x - y) for x, y in zip(ds(a[k]), ds(b[k]))) for k in diff if ds(a[k]) and len(ds(a[k])) == len(ds(b[k]))), default=0)
print(f"GPU vs CPU: {len(shared)} variants scored by both, {len(shared) - len(diff)} identical, {len(diff)} differ "
      f"(largest delta-score difference {dmax:.2f})")
PY
fi
