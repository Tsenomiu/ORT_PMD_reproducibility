#!/usr/bin/env bash
# Run four datasets across 22 one-chromosome GLIMPSE2 comparisons.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
manifest=${1:?Usage: run_per_chromosome_concordance.sh MANIFEST OUTPUT_DIR}
outdir=${2:?}
mkdir -p "$outdir"

while IFS=$'\t' read -r individual treatment comparator target; do
  [[ -z "$individual" || "$individual" == individual || "$individual" == \#* ]] && continue
  for chrom in $(seq 1 22); do
    freq=${CONCORDANCE_FREQ_PATTERN//\{CHROM\}/$chrom}
    for path in "$freq" "$comparator" "$target"; do
      case "$path" in "$CONCORDANCE_INPUT_ROOT"/*) ;; *) echo "$path is outside CONCORDANCE_INPUT_ROOT" >&2; exit 2;; esac
    done
    list="$outdir/list_${individual}_${treatment}_chr${chrom}.txt"
    printf '%s\t/data/%s\t/data/%s\t/data/%s\n' "$chrom" \
      "${freq#"$CONCORDANCE_INPUT_ROOT"/}" "${comparator#"$CONCORDANCE_INPUT_ROOT"/}" \
      "${target#"$CONCORDANCE_INPUT_ROOT"/}" > "$list"
    "$DOCKER" run --rm -v "$CONCORDANCE_INPUT_ROOT:/data:ro" -v "$outdir:/out" \
      "$GLIMPSE2_IMAGE" GLIMPSE2_concordance \
      --input "/out/$(basename "$list")" --gt-val --af-tag AF \
      --bins 0 0.001 0.01 0.02 0.05 0.10 0.20 0.30 0.40 0.50 \
      --thread 6 --output "/out/pc_${individual}_${treatment}_chr${chrom}"
  done
done < "$manifest"
"$script_dir/harvest_perchrom_gcsv.sh" "$manifest" "$outdir" "$outdir/perchr_gcsv_rows.txt"
