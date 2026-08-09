#!/usr/bin/env bash
# Merge sequencing units and remove PCR duplicates per read set.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

manifest=${1:?Usage: merge_markdup.sh MANIFEST OUTPUT_DIR}
outdir=${2:?}
mkdir -p "$outdir"

output_ids=()
while IFS= read -r output_id; do
  output_ids[${#output_ids[@]}]="$output_id"
done < <(awk -F '\t' 'NR>1 && $1!=""{print $1}' "$manifest" | sort -u)
for output_id in "${output_ids[@]}"; do
  mode=$(awk -F '\t' -v id="$output_id" 'NR>1 && $1==id{print $2}' "$manifest" | sort -u)
  [[ "$mode" == SE || "$mode" == PE ]] || { echo "One SE/PE mode required for $output_id" >&2; exit 1; }
  list="$outdir/${output_id}.bam.list"
  awk -F '\t' -v id="$output_id" 'NR>1 && $1==id{print $3}' "$manifest" > "$list"
  while IFS= read -r bam; do [[ -s "$bam" ]] || { echo "Missing $bam" >&2; exit 1; }; done < "$list"
  merged="$outdir/${output_id}.merged.bam"
  "$SAMTOOLS" merge -f -b "$list" "$merged"
  if [[ "$mode" == PE ]]; then
    "$SAMTOOLS" sort -n -@ "$THREADS" -o "$outdir/${output_id}.name.bam" "$merged"
    "$SAMTOOLS" fixmate -m "$outdir/${output_id}.name.bam" "$outdir/${output_id}.fixmate.bam"
    "$SAMTOOLS" sort -@ "$THREADS" -o "$outdir/${output_id}.positions.bam" "$outdir/${output_id}.fixmate.bam"
    "$SAMTOOLS" markdup -@ "$THREADS" -r "$outdir/${output_id}.positions.bam" "$outdir/${output_id}.dedup.bam"
    rm -f "$outdir/${output_id}.name.bam" "$outdir/${output_id}.fixmate.bam"
  else
    "$SAMTOOLS" sort -@ "$THREADS" -o "$outdir/${output_id}.positions.bam" "$merged"
    "$SAMTOOLS" markdup -s -@ "$THREADS" -r "$outdir/${output_id}.positions.bam" "$outdir/${output_id}.dedup.bam"
  fi
  "$SAMTOOLS" index "$outdir/${output_id}.dedup.bam"
  rm -f "$merged" "$outdir/${output_id}.positions.bam" "$list"
done
