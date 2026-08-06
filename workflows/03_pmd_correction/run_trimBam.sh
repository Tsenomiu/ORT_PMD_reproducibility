#!/usr/bin/env bash
# Soft-clip both read ends while retaining clipped sequence in the BAM.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
mode=${1:?Usage: run_trimBam.sh SE|PE LENGTH INPUT_BAM OUTPUT_BAM}
length=${2:?}
input=${3:?}
output=${4:?}
mkdir -p "$(dirname "$output")"
if [[ "$mode" == SE ]]; then
  "$TRIMBAM" trimBam "$input" - "$length" -c | "$SAMTOOLS" sort -@ "$THREADS" -o "$output" -
elif [[ "$mode" == PE ]]; then
  "$TRIMBAM" trimBam "$input" - "$length" -c \
    | "$SAMTOOLS" sort -n -@ "$THREADS" -o "$output.name.bam" -
  "$SAMTOOLS" fixmate -m "$output.name.bam" "$output.fixmate.bam"
  "$SAMTOOLS" sort -@ "$THREADS" -o "$output" "$output.fixmate.bam"
  rm -f "$output.name.bam" "$output.fixmate.bam"
else
  echo "MODE must be SE or PE" >&2; exit 2
fi
"$SAMTOOLS" index "$output"
