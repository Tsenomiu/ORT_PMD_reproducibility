#!/usr/bin/env bash
# Publication wrapper around the exact recovered primary downstream commands.
# The original BAM-to-genotype command/random state is not claimed or recreated.
set -euo pipefail

input_prefix=${1:?Usage: run_primary_from_processed_plink.sh INPUT_PLINK_PREFIX OUTPUT_DIR}
output_dir=${2:?Usage: run_primary_from_processed_plink.sh INPUT_PLINK_PREFIX OUTPUT_DIR}

: "${READ2_PY:?Set READ2_PY to READv2 tag 2.1.0 READ2.py}"
PLINK=${PLINK:-plink}
PYTHON=${PYTHON:-python3}

for suffix in bed bim fam; do
  [[ -s "${input_prefix}.${suffix}" ]] || {
    echo "Missing processed primary input: ${input_prefix}.${suffix}" >&2
    exit 1
  }
done

[[ $(wc -l < "${input_prefix}.fam") -eq 27 ]] || {
  echo "The preserved primary boundary must contain ORT1–ORT27 (27 people)." >&2
  exit 1
}

read2_dir=$(cd "$(dirname "$READ2_PY")" && pwd -P)
READ2_PY="${read2_dir}/$(basename "$READ2_PY")"
[[ -s "$READ2_PY" ]] || { echo "READ2.py not found: $READ2_PY" >&2; exit 1; }

mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd -P)
for suffix in bed bim fam; do
  [[ ! -e "$output_dir/filtered.${suffix}" ]] || {
    echo "Refusing to overwrite $output_dir/filtered.${suffix}" >&2
    exit 1
  }
done

# Exact recovered filtering options.
"$PLINK" --bfile "$input_prefix" --mind 0.8 --make-bed --out "$output_dir/filtered"

[[ $(wc -l < "$output_dir/filtered.fam") -eq 21 ]] || {
  echo "Expected 21 people after --mind 0.8." >&2
  exit 1
}

removed=$(
  comm -23 \
    <(awk '{print $2}' "${input_prefix}.fam" | LC_ALL=C sort) \
    <(awk '{print $2}' "$output_dir/filtered.fam" | LC_ALL=C sort) \
    | paste -sd, -
)
[[ "$removed" == "ORT1,ORT17,ORT19,ORT20,ORT22,ORT25" ]] || {
  echo "Unexpected --mind exclusions: $removed" >&2
  exit 1
}

# Exact recovered READv2 options; default median normalization is retained.
(
  cd "$output_dir"
  "$PYTHON" "$READ2_PY" -i filtered --window_size 500000
)

echo "Primary READv2 run complete: $output_dir"
