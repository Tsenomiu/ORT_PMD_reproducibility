#!/usr/bin/env bash
# Run the fixed-comparator 22-region GLIMPSE2 design.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

individual=${1:?Usage: run_fixed_comparator_concordance.sh IND TREAT COMPARATOR_VCF TARGET_VCF OUTPUT_DIR}
treatment=${2:?}
comparator=${3:?}
target=${4:?}
outdir=${5:?}
mkdir -p "$outdir"

case "$comparator" in "$CONCORDANCE_INPUT_ROOT"/*) ;; *) echo "Comparator must be below CONCORDANCE_INPUT_ROOT" >&2; exit 2;; esac
withds="$outdir/${individual}_${treatment}.withDS.vcf.gz"
"$PYTHON" "$script_dir/add_ds_from_gp.py" "$target" | "$BGZIP" -c > "$withds"
"$TABIX" -f -p vcf "$withds"

host_list="$outdir/${individual}_${treatment}.regions.container.txt"
container_list="/out/$(basename "$host_list")"
: > "$host_list"
for chrom in $(seq 1 22); do
  freq=${CONCORDANCE_FREQ_PATTERN//\{CHROM\}/$chrom}
  case "$freq" in "$CONCORDANCE_INPUT_ROOT"/*) ;; *) echo "Frequency VCF must be below CONCORDANCE_INPUT_ROOT" >&2; exit 2;; esac
  printf '%s\t/data/%s\t/data/%s\t/out/%s\n' \
    "$chrom" "${freq#"$CONCORDANCE_INPUT_ROOT"/}" \
    "${comparator#"$CONCORDANCE_INPUT_ROOT"/}" "$(basename "$withds")" >> "$host_list"
done

"$DOCKER" run --rm \
  -v "$CONCORDANCE_INPUT_ROOT:/data:ro" -v "$outdir:/out" \
  "$GLIMPSE2_IMAGE" GLIMPSE2_concordance \
  --input "$container_list" --gt-val --af-tag AF \
  --bins 0 0.001 0.01 0.02 0.05 0.10 0.20 0.30 0.40 0.50 \
  --thread 16 --output "/out/conc_${individual}_${treatment}"
