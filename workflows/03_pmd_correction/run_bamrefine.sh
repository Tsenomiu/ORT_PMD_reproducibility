#!/usr/bin/env bash
# Run bamRefine with a configured SNP list and terminal threshold.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
length=${1:?Usage: run_bamrefine.sh 5|10 INPUT_BAM OUTPUT_BAM}
input=${2:?}
output=${3:?}
[[ "$length" == 5 || "$length" == 10 ]] || { echo "Threshold must be 5 or 10" >&2; exit 2; }
mkdir -p "$(dirname "$output")"
"$BAMREFINE" --snps "$BAMREFINE_SNP_LIST" --pmd-length-threshold "$length" \
  --threads "$THREADS" "$input" "$output"
