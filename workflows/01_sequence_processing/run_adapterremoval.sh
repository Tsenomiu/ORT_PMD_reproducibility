#!/usr/bin/env bash
# Run AdapterRemoval v2.2.2 with the study parameters.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

manifest=${1:?Usage: run_adapterremoval.sh FASTQ_MANIFEST OUTPUT_DIR}
outdir=${2:?Usage: run_adapterremoval.sh FASTQ_MANIFEST OUTPUT_DIR}
mkdir -p "$outdir"

while IFS=$'\t' read -r sample library lane r1 r2; do
  [[ -z "$sample" || "$sample" == sample || "$sample" == \#* ]] && continue
  [[ -s "$r1" && -s "$r2" ]] || { echo "Missing FASTQ for $sample/$library/$lane" >&2; exit 1; }
  prefix="$outdir/${sample}_${library}_${lane}"
  "$ADAPTERREMOVAL" \
    --file1 "$r1" --file2 "$r2" \
    --basename "$prefix" \
    --trimns --trimqualities \
    --minlength 30 --minquality 25 \
    --collapse --gzip
done < "$manifest"
