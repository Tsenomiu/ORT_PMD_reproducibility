#!/bin/sh
set -eu

CODE_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
RELEASE_ROOT=$(CDPATH= cd -- "$CODE_DIR/.." && pwd)
BUILD_DIR=${1:-"$CODE_DIR/reproduced_figures"}
PYTHON=${PYTHON:-python3}
SHARED="$CODE_DIR/_shared/ancibd_v62"
SUMMARY="$RELEASE_ROOT/data/summary"
MPLCONFIGDIR=${MPLCONFIGDIR:-"$BUILD_DIR/.mplconfig"}
export MPLCONFIGDIR

mkdir -p "$BUILD_DIR" "$MPLCONFIGDIR"
"$PYTHON" -c 'import matplotlib, numpy, PIL'

check_pdf() {
  label=$1
  path=$2
  test -s "$path"
  if command -v pdfinfo >/dev/null 2>&1; then
    pages=$(pdfinfo "$path" | awk '/^Pages:/ {print $2}')
    test "$pages" = 1
  fi
  printf 'PASS: %s -> %s\n' "$label" "$path"
}

if "$CODE_DIR/figure_01_workflow/render_svg.sh" "$BUILD_DIR/figure_01_workflow.pdf"; then
  check_pdf 'Figure 1' "$BUILD_DIR/figure_01_workflow.pdf"
else
  status=$?
  if test "$status" -eq 77; then
    printf '%s\n' 'SKIP: Figure 1 (no supported SVG-to-PDF converter)'
  else
    exit "$status"
  fi
fi

"$PYTHON" "$CODE_DIR/figure_02_alt_fraction/make_figure.py" \
  "$SUMMARY/damage/alt_fraction_recomputed_dp3_adsum.csv" \
  "$BUILD_DIR/figure_02_alt_fraction.pdf"
check_pdf 'Figure 2' "$BUILD_DIR/figure_02_alt_fraction.pdf"

"$PYTHON" "$CODE_DIR/figure_03_imputation/make_figure.py" \
  --summary "$SUMMARY/imputation/imputation_summary.csv" \
  --output "$BUILD_DIR/figure_03_imputation.pdf"
check_pdf 'Figure 3' "$BUILD_DIR/figure_03_imputation.pdf"

"$PYTHON" "$CODE_DIR/figure_04_ibd_summary/make_figures.py" \
  --asymmetric "$SUMMARY/ancibd/asymmetric_fulludg_raw.tsv" \
  --summary "$SUMMARY/ancibd/cross_treatment_pair_summary.tsv" \
  --resource-manifest "$SHARED/resource_manifest_minimal.json" \
  --output-dir "$BUILD_DIR"
check_pdf 'Figure 4' "$BUILD_DIR/figure_04_ibd_summary.pdf"
check_pdf 'Supplementary Figure S6' "$BUILD_DIR/supplementary_06_ibd_length_count.pdf"

"$PYTHON" "$CODE_DIR/figure_05_ibd_karyogram/make_figure.py" \
  --base "$SHARED" \
  --resource-manifest "$SHARED/resource_manifest_minimal.json" \
  --output-dir "$BUILD_DIR"
check_pdf 'Figure 5' "$BUILD_DIR/figure_05_ibd_karyogram.pdf"

"$PYTHON" "$CODE_DIR/figure_06_kinship/make_figure.py" \
  --output "$BUILD_DIR/figure_06_kinship.pdf"
check_pdf 'Figure 6' "$BUILD_DIR/figure_06_kinship.pdf"

ORT_FIGURE_OUTPUT="$BUILD_DIR/figure_07_pca.pdf" \
  "$PYTHON" "$CODE_DIR/figure_07_pca/make_figure.py"
check_pdf 'Figure 7 (aggregate reference background)' "$BUILD_DIR/figure_07_pca.pdf"

"$PYTHON" "$CODE_DIR/supplementary_01_pmd_uncorrected/make_figure.py" \
  --output "$BUILD_DIR/supplementary_01_pmd_uncorrected.pdf"
check_pdf 'Supplementary Figure S1' "$BUILD_DIR/supplementary_01_pmd_uncorrected.pdf"

printf '%s\n' 'SKIP: Supplementary Figure S2 (controlled archaeological coordinates and contour source)'
printf '%s\n' 'SKIP: Supplementary Figure S3 (controlled VOX volumes and selected slice coordinates)'

"$PYTHON" "$CODE_DIR/supplementary_04_pmd_corrected/make_figure.py" \
  --output "$BUILD_DIR/supplementary_04_pmd_corrected.pdf"
check_pdf 'Supplementary Figure S4' "$BUILD_DIR/supplementary_04_pmd_corrected.pdf"

"$PYTHON" "$CODE_DIR/supplementary_05_ibd_matrix/make_figure.py" \
  --summary "$SUMMARY/ancibd/cross_treatment_pair_summary.tsv" \
  --resource-manifest "$SHARED/resource_manifest_minimal.json" \
  --output "$BUILD_DIR/supplementary_05_ibd_matrix.pdf"
check_pdf 'Supplementary Figure S5' "$BUILD_DIR/supplementary_05_ibd_matrix.pdf"

if command -v pdftex >/dev/null 2>&1; then
  ORT_FIGURE_OUTPUT_DIR="$BUILD_DIR" \
    "$PYTHON" "$CODE_DIR/supplementary_07_08_pca_titration/make_figures.py"
  check_pdf 'Supplementary Figure S7 (aggregate reference background)' "$BUILD_DIR/supplementary_07_pca_titration_ort15.pdf"
  check_pdf 'Supplementary Figure S8 (aggregate reference background)' "$BUILD_DIR/supplementary_08_pca_titration_ort16.pdf"
else
  printf '%s\n' 'SKIP: Supplementary Figures S7-S8 (pdftex not installed)'
fi

printf '%s\n' 'PASS: all figure generators supported by the included inputs and installed tools completed.'
printf '%s\n' 'NOTE: Figure 7 and S7-S8 use aggregate reference backgrounds; exact individual-level backgrounds require external AADR data.'
printf 'Rebuilt files: %s\n' "$BUILD_DIR"
