#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Reconstructed from successful project logs. This is not claimed as the
# historical top-level script. It reruns the deterministic plink2tkrelated
# stage from the retained, checksummed PED/MAP inputs.
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_tkgwv2_from_processed_peds.sh \
    --tkgwv2 TKGWV2.py \
    --frequency-file EAS_FREQUENCY_FILE \
    --ped-map-dir PED_MAP_DIRECTORY \
    --output-dir NEW_OUTPUT_DIRECTORY \
    [--python-bin python3]
USAGE
}

tkgwv2=""
frequency_file=""
ped_map_dir=""
output_dir=""
python_bin="${PYTHON:-python3}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tkgwv2) tkgwv2=${2:-}; shift 2 ;;
    --frequency-file) frequency_file=${2:-}; shift 2 ;;
    --ped-map-dir) ped_map_dir=${2:-}; shift 2 ;;
    --output-dir) output_dir=${2:-}; shift 2 ;;
    --python-bin) python_bin=${2:-}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

for value in "$tkgwv2" "$frequency_file" "$ped_map_dir" "$output_dir"; do
  [[ -n "$value" ]] || { usage; exit 2; }
done

[[ -f "$tkgwv2" ]] || { echo "TKGWV2.py not found: $tkgwv2" >&2; exit 3; }
[[ -f "$frequency_file" ]] || { echo "Frequency file not found" >&2; exit 3; }
[[ -d "$ped_map_dir" ]] || { echo "PED/MAP directory not found" >&2; exit 3; }
command -v "$python_bin" >/dev/null 2>&1 || {
  echo "Python command not found: $python_bin" >&2
  exit 3
}
[[ ! -e "$output_dir" ]] || {
  echo "Refusing to overwrite existing output: $output_dir" >&2
  exit 4
}

tkgwv2=$(cd "$(dirname "$tkgwv2")" && pwd -P)/$(basename "$tkgwv2")
frequency_file=$(cd "$(dirname "$frequency_file")" && pwd -P)/$(basename "$frequency_file")
ped_map_dir=$(cd "$ped_map_dir" && pwd -P)
package_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)

fu=(
  ORT15_fu_raw ORT15_fu_bamrefine5 ORT15_fu_rescale5 ORT15_fu_trim5
  ORT16_fu_raw ORT16_fu_bamrefine5 ORT16_fu_rescale5 ORT16_fu_trim5
)
nu=(
  ORT15_nu_raw ORT15_nu_bamrefine5 ORT15_nu_rescale5 ORT15_nu_trim5
  ORT16_nu_raw ORT16_nu_bamrefine5 ORT16_nu_rescale5 ORT16_nu_trim5
)

for sample in "${fu[@]}" "${nu[@]}"; do
  for extension in ped map; do
    [[ -f "$ped_map_dir/$sample.$extension" ]] || {
      echo "Required input not found: $sample.$extension" >&2
      exit 3
    }
  done
done

mkdir -p "$output_dir"/{fu_fu,nu_nu,nu_fu}
output_dir=$(cd "$output_dir" && pwd -P)

for sample in "${fu[@]}"; do
  for run in fu_fu nu_fu; do
    ln -s "$ped_map_dir/$sample.ped" "$output_dir/$run/$sample.ped"
    ln -s "$ped_map_dir/$sample.map" "$output_dir/$run/$sample.map"
  done
done
for sample in "${nu[@]}"; do
  for run in nu_nu nu_fu; do
    ln -s "$ped_map_dir/$sample.ped" "$output_dir/$run/$sample.ped"
    ln -s "$ped_map_dir/$sample.map" "$output_dir/$run/$sample.map"
  done
done

for run in fu_fu nu_nu nu_fu; do
  cp "$package_dir/dyads/$run.tsv" "$output_dir/$run/dyads.txt"
  (
    cd "$output_dir/$run"
    "$python_bin" "$tkgwv2" plink2tkrelated \
      --freqFile "$frequency_file" \
      --dyads dyads.txt > run.stdout 2>&1
  )
done

"$python_bin" "$package_dir/validation/validate_results.py" \
  --results-dir "$output_dir" \
  --output "$output_dir/validation_results.tsv"

echo "Completed and validated the deterministic processed-PED rerun: $output_dir"
