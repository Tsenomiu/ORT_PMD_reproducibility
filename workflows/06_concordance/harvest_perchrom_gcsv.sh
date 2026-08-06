#!/usr/bin/env bash
# Collect one GCsV record from each per-chromosome output.
set -euo pipefail
manifest=${1:?Usage: harvest_perchrom_gcsv.sh MANIFEST CONCORDANCE_DIR OUTPUT}
indir=${2:?}
output=${3:?}
: > "$output"
while IFS=$'\t' read -r individual treatment _ _; do
  [[ -z "$individual" || "$individual" == "individual" || "$individual" == \#* ]] && continue
  for chrom in $(seq 1 22); do
    file="$indir/pc_${individual}_${treatment}_chr${chrom}.error.spl.txt.gz"
    [[ -s "$file" ]] || { echo "Missing $file" >&2; exit 1; }
    row=$(gzip -dc "$file" | awk '$1=="GCsV"{print; exit}')
    [[ -n "$row" ]] || { echo "No GCsV row in $file" >&2; exit 1; }
    printf '%s %s %s | %s\n' "$individual" "$treatment" "$chrom" "$row" >> "$output"
  done
done < "$manifest"
[[ $(wc -l < "$output") -eq 88 ]] || { echo "Expected 88 GCsV rows" >&2; exit 1; }
