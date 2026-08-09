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

printf '%s\n' '== READv2 ORT15–ORT16 consistency =='
"$PYTHON" - "$ROOT/data/summary/kinship/readv2_ort15_ort16.tsv" \
  "$ROOT/figures/figure_06_kinship/input_values.tsv" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = {row["analysis"]: row for row in csv.DictReader(handle, delimiter="\t")}
with open(sys.argv[2], newline="", encoding="utf-8") as handle:
    values = {r["metric"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}

required = {"primary_chr1_22_X_Y", "autosome_only_sensitivity"}
if set(rows) != required:
    raise SystemExit(f"Unexpected READv2 analyses: {sorted(rows)}")

primary = rows["primary_chr1_22_X_Y"]
sensitivity = rows["autosome_only_sensitivity"]
observed = float(primary["KinshipCoefficient"])
plotted = float(values["READv2_normalized_kinship"])
if not math.isclose(observed, plotted, rel_tol=0.0, abs_tol=1e-15):
    raise SystemExit(f"READv2 KC mismatch: summary={observed}, figure={plotted}")
if int(primary["OverlapNSNPs"]) != 169172 or primary["First_degree_subtype"] != "Parent-offspring":
    raise SystemExit("Primary READv2 row does not match the archived chr1-22/X/Y result")
if (
    int(sensitivity["OverlapNSNPs"]) != 157451
    or sensitivity["First_degree_subtype"] != "Parent-offspring"
    or not math.isclose(float(sensitivity["KinshipCoefficient"]), 0.23159213948114488,
                        rel_tol=0.0, abs_tol=1e-15)
):
    raise SystemExit("Autosome-only READv2 sensitivity row is inconsistent")
PY
printf '%s\n' 'PASS  READv2 primary and autosome-only sensitivity rows are consistent'

printf '%s\n' '== Reported-result summary validation =='
"$PYTHON" "$ROOT/validate_reported_results.py"

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

printf '%s\n' '== Figure regeneration =='
PYTHON="$PYTHON" bash "$ROOT/figures/reproduce_figures.sh" "$BUILD_DIR/figures"

printf '\n%s\n' 'Reproduction complete: all locally reproducible checks passed.'
printf 'Build products: %s\n' "$BUILD_DIR"
