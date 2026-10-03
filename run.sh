#!/usr/bin/env bash
# Run brava-annotate with whatever is available: Apptainer/Singularity, Docker, or pixi.
#
#   ./run.sh setup -r brava_resources                 # once (~5 GB)
#   ./run.sh -r brava_resources -t 16 -o out chr*.vcf.gz
#
# Env: BRAVA_IMAGE (container image), BRAVA_SIF (local .sif), BRAVA_RUNNER (apptainer|docker|pixi)
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
IMAGE=${BRAVA_IMAGE:-ghcr.io/brava-genetics/variant-annotation:latest}

gpu=0
for a in "$@"; do [[ $a == --gpu ]] && gpu=1; done
(( gpu )) && IMAGE=${BRAVA_IMAGE:-${IMAGE%:*}:latest-gpu}

# Directories the container must see: cwd, and the parent of every path-like argument
binds=("$PWD")
prev=""
for a in "$@"; do
  if [[ $prev == -o || $prev == --out || $prev == -r || $prev == --resources ]]; then mkdir -p "$a"; fi
  if [[ -e $a ]]; then binds+=("$(cd "$(dirname "$a")" && pwd)"); fi
  prev=$a
done
[[ -n ${BRAVA_RESOURCES:-} ]] && { mkdir -p "$BRAVA_RESOURCES"; binds+=("$(cd "$BRAVA_RESOURCES" && pwd)"); }
mapfile -t binds < <(printf '%s\n' "${binds[@]}" | sort -u)

runner=${BRAVA_RUNNER:-}
if [[ -z $runner ]]; then
  if command -v vep >/dev/null && command -v spliceai >/dev/null; then runner=native
  elif command -v apptainer >/dev/null || command -v singularity >/dev/null; then runner=apptainer
  elif command -v docker >/dev/null; then runner=docker
  elif command -v pixi >/dev/null; then runner=pixi
  else echo "Need one of: apptainer/singularity, docker, or pixi (https://pixi.sh)" >&2; exit 1; fi
fi

case $runner in
  native) exec "$HERE/bin/brava-annotate" "$@" ;;
  apptainer)
    exe=$(command -v apptainer || command -v singularity)
    flags=(--cleanenv --env "BRAVA_RESOURCES=${BRAVA_RESOURCES:-}")
    (( gpu )) && flags+=(--nv)
    exec "$exe" exec "${flags[@]}" --bind "$(IFS=,; echo "${binds[*]}")" \
      "${BRAVA_SIF:-docker://$IMAGE}" brava-annotate "$@" ;;
  docker)
    flags=(--rm -u "$(id -u):$(id -g)" -w "$PWD" -e "BRAVA_RESOURCES=${BRAVA_RESOURCES:-}")
    (( gpu )) && flags+=(--gpus all)
    for b in "${binds[@]}"; do flags+=(-v "$b:$b"); done
    exec docker run "${flags[@]}" "$IMAGE" brava-annotate "$@" ;;
  pixi)
    env=default; (( gpu )) && env=gpu
    exec pixi run --manifest-path "$HERE/pixi.toml" -e "$env" "$HERE/bin/brava-annotate" "$@" ;;
  *) echo "unknown BRAVA_RUNNER=$runner" >&2; exit 1 ;;
esac
