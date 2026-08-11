#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BUILD_DIR=${ORT_BUILD_DIR:-"$ROOT/_build"}
PYTHON=${PYTHON:-python3}
export PYTHONDONTWRITEBYTECODE=1

mkdir -p "$BUILD_DIR/tables" "$BUILD_DIR/figures"

printf '%s\n' '== Repository checks =='
bash "$ROOT/check_repository.sh"

printf '%s\n' '== Concordance column-identity unit test =='
"$PYTHON" "$ROOT/workflows/06_concordance/test_summarize_concordance.py"
printf '%s\n' 'PASS  concordance dosage/best-guess column identity'

printf '%s\n' '== Posterior-dosage reconstruction unit test =='
"$PYTHON" "$ROOT/workflows/06_concordance/test_add_ds_from_gp.py"
printf '%s\n' 'PASS  dosage reconstruction and four-decimal precision'

printf '%s\n' '== Chromosome-jackknife reconstruction =='
JACK_DIR="$BUILD_DIR/tables/jackknife"
"$PYTHON" "$ROOT/workflows/06_concordance/compute_jackknife.py" \
  "$ROOT/data/summary/imputation/perchr_gcsv_rows.txt" \
  "$JACK_DIR"
cmp "$ROOT/data/summary/imputation/perchr_nrd_counts.csv" \
  "$JACK_DIR/perchr_nrd_counts.csv"
cmp "$ROOT/data/summary/imputation/jackknife_deltaNRD.csv" \
  "$JACK_DIR/jackknife_deltaNRD.csv"
printf '%s\n' 'PASS  jackknife outputs reproduce byte-for-byte'

printf '%s\n' '== Reported-result summary validation =='
"$PYTHON" "$ROOT/validate_reported_results.py"

printf '%s\n' '== Recovered TKGWV2 validation =='
"$PYTHON" "$ROOT/workflows/08_kinship/tkgwv2/validation/validate_results.py"

printf '%s\n' '== Recovered READv2 validation =='
"$PYTHON" "$ROOT/workflows/08_kinship/readv2/validate_reported_results.py" \
  --output "$BUILD_DIR/tables/readv2_validation_results.tsv"

printf '%s\n' '== Recovered Table 1 validation =='
"$PYTHON" "$ROOT/workflows/06_concordance/table1/validate_table1.py"

printf '%s\n' '== Main-output source coverage =='
"$PYTHON" "$ROOT/docs/validation/validate_main_output_sources.py"

printf '%s\n' '== Public PCA-coordinate verification =='
"$PYTHON" "$ROOT/workflows/10_pca/verify_pca.py" \
  "$ROOT/data/summary/pca/pca_ort_queries_matched48.evec" \
  "$ROOT/data/summary/pca/pca_aadr_matched48.eval" \
  "$ROOT/data/summary/pca/query_manifest_matched48.tsv" \
  "$ROOT/data/summary/pca/reference_population_aggregates.tsv" \
  "$BUILD_DIR/tables/pca_metrics_verified.csv" \
  --reference "$ROOT/data/summary/pca/pca_metrics_matched48.csv" \
  --absolute-tolerance 1e-12
printf '%s\n' 'PASS  public PCA metrics reproduce within 1e-12 absolute tolerance'

printf '%s\n' '== Recovered PCA workflow and raster validation =='
if command -v pdftex >/dev/null 2>&1 && command -v pdftoppm >/dev/null 2>&1; then
  PYTHON="$PYTHON" bash "$ROOT/workflows/10_pca/recovery/run_validation_and_figures.sh" \
    "$BUILD_DIR/pca_recovery"
else
  printf '%s\n' 'SKIP  recovered PCA raster validation (pdftex/pdftoppm unavailable)'
fi

printf '%s\n' '== Figure regeneration =='
PYTHON="$PYTHON" bash "$ROOT/figures/reproduce_figures.sh" "$BUILD_DIR/figures"

printf '\n%s\n' 'Reproduction complete: all locally reproducible checks passed.'
printf 'Build products: %s\n' "$BUILD_DIR"
