#!/usr/bin/env bash
# Run mapDamage2 rescaling; LENGTH is 5, 10, or default12.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
length=${1:?Usage: run_mapdamage.sh 5|10|default12 INPUT_BAM OUTPUT_DIR}
input=${2:?}
outdir=${3:?}
mkdir -p "$outdir"
args=(-i "$input" -r "$HS37D5_FASTA" -d "$outdir" --rescale --merge-libraries)
if [[ "$length" != default12 ]]; then args+=(--seq-length "$length"); fi
"$MAPDAMAGE" "${args[@]}"
