#!/usr/bin/env bash
# Recalculate PCA metrics and reproduce the three manuscript PCA rasters.
set -euo pipefail
export LC_ALL=C

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
repo_root=$(cd "$script_dir/../../.." && pwd -P)
build_dir=${1:-"$repo_root/_build/pca_recovery"}
PYTHON=${PYTHON:-python3}
PDFTOPPM=${PDFTOPPM:-pdftoppm}

mkdir -p "$build_dir"
build_dir=$(cd "$build_dir" && pwd -P)
export MPLCONFIGDIR=${MPLCONFIGDIR:-"$build_dir/.mplconfig"}
mkdir -p "$MPLCONFIGDIR"

command -v "$PDFTOPPM" >/dev/null 2>&1 || {
  echo "pdftoppm is required for the 150-dpi raster comparison" >&2
  exit 2
}
command -v pdftex >/dev/null 2>&1 || {
  echo "pdftex is required to rebuild Supplementary Figures S7 and S8" >&2
  exit 2
}
"$PYTHON" -c 'import matplotlib, numpy, PIL'

ORT_FIGURE_OUTPUT="$build_dir/figure_07_pca.pdf" \
  "$PYTHON" "$repo_root/figures/figure_07_pca/make_figure.py"
ORT_FIGURE_OUTPUT_DIR="$build_dir" \
  "$PYTHON" "$repo_root/figures/supplementary_07_08_pca_titration/make_figures.py"

"$PDFTOPPM" -r 150 -png -singlefile "$build_dir/figure_07_pca.pdf" \
  "$build_dir/main_pca_raster_150dpi"
"$PDFTOPPM" -r 150 -png -singlefile \
  "$build_dir/supplementary_07_pca_titration_ort15.pdf" \
  "$build_dir/titration_ort15_raster_150dpi"
"$PDFTOPPM" -r 150 -png -singlefile \
  "$build_dir/supplementary_08_pca_titration_ort16.pdf" \
  "$build_dir/titration_ort16_raster_150dpi"

"$PYTHON" "$script_dir/validate_pca_recovery.py" \
  --raster-dir "$build_dir" \
  --output "$build_dir/validation_results_recomputed.tsv" \
  --absolute-tolerance 1e-12

echo "PASS: three PCA PDFs rebuilt and their 150-dpi rasters matched"
echo "Validation outputs: $build_dir"
