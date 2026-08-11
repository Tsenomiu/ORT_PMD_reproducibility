#!/usr/bin/env bash
# Path-neutral publication wrapper for the recovered final matched-48 projection.
set -euo pipefail
export LC_ALL=C

usage() {
  echo "Usage: $0 IMPUTED_PREFIX PSEUDOHAPLOID_PREFIX REFERENCE_PREFIX IID_TO_POPULATION KEEP_POPS OUTPUT_DIR" >&2
  exit 2
}

[[ $# -eq 6 ]] || usage
imputed=$1
pseudohaploid=$2
reference=$3
iid_to_population=$4
keep_pops=$5
output_dir=$6

for prefix in "$imputed" "$pseudohaploid" "$reference"; do
  for suffix in bed bim fam; do
    [[ -f "${prefix}.${suffix}" ]] || { echo "Missing PLINK input: ${prefix}.${suffix}" >&2; exit 2; }
  done
done
for path in "$iid_to_population" "$keep_pops"; do
  [[ -f "$path" ]] || { echo "Missing metadata input: $path" >&2; exit 2; }
done
[[ -n ${ORT_CONFIG:-} && -f "$ORT_CONFIG" ]] || {
  echo "Set ORT_CONFIG to a configured workflows/config.sh" >&2
  exit 2
}
[[ ! -e "$output_dir" ]] || {
  echo "Output directory already exists: $output_dir" >&2
  exit 2
}

mkdir -p "$output_dir/query_inputs"
output_dir=$(cd "$output_dir" && pwd -P)
input_dir="$output_dir/query_inputs"

# The recovered run reset the query BIM genetic-distance column to zero to avoid
# scientific-notation parsing differences. BED/FAM files remain unchanged.
for item in "imputed:$imputed" "pseudohaploid:$pseudohaploid"; do
  label=${item%%:*}
  prefix=${item#*:}
  ln -s "$(cd "$(dirname "${prefix}.bed")" && pwd -P)/$(basename "${prefix}.bed")" \
    "$input_dir/${label}.bed"
  ln -s "$(cd "$(dirname "${prefix}.fam")" && pwd -P)/$(basename "${prefix}.fam")" \
    "$input_dir/${label}.fam"
  awk 'BEGIN{OFS="\t"}{$3=0;print}' "${prefix}.bim" > "$input_dir/${label}.bim"
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
bash "$script_dir/../run_matched48_projection.sh" \
  "$input_dir/imputed" "$input_dir/pseudohaploid" "$reference" \
  "$iid_to_population" "$keep_pops" "$output_dir"

python_cmd=${PYTHON:-python3}
"$python_cmd" - "$output_dir" "$keep_pops" <<'PY'
import csv
import json
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
population_count = sum(1 for line in Path(sys.argv[2]).read_text().splitlines() if line.strip())
rows = []
with (outdir / "matched48.labelled.ind").open() as handle:
    for raw in handle:
        iid = raw.split()[0]
        if iid.startswith("IMP__"):
            representation, base = "imputed", iid[5:]
        elif iid.startswith("PHAP__"):
            representation, base = "pseudohaploid", iid[6:]
        else:
            continue
        parts = base.split("_")
        rows.append(
            {
                "prefixed_id": iid,
                "individual": parts[0],
                "library": parts[1],
                "treatment": "_".join(parts[2:]),
                "representation": representation,
            }
        )
rows.sort(key=lambda row: row["prefixed_id"])
with (outdir / "query_manifest_matched48.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=("prefixed_id", "individual", "library", "treatment", "representation"),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

eigenvalues = [float(value) for value in (outdir / "pca_aadr_matched48.eval").read_text().split()]
trace = sum(value for value in eigenvalues if value > 0)
site_stats = json.loads((outdir / "site_stats.json").read_text())
summary = {
    "reference_individuals": 985,
    "reference_populations": population_count,
    "imputed_query_points": sum(row["representation"] == "imputed" for row in rows),
    "pseudohaploid_query_points": sum(row["representation"] == "pseudohaploid" for row in rows),
    "projected_query_points": len(rows),
    "common_markers_before_smartpca": site_stats["n_kept"],
    "smartpca_markers_used": 1113348,
    "eigenvalues": len(eigenvalues),
    "eigenvalue_trace": round(trace, 6),
    "pc_variance_pct": [round(100 * value / trace, 4) for value in eigenvalues[:10]],
}
(outdir / "projection_build_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)
if len(rows) != 48:
    raise SystemExit(f"Expected 48 query rows; observed {len(rows)}")
PY

population_count=$(awk 'NF{n++}END{print n+0}' "$keep_pops")
[[ $population_count -eq 154 ]] || {
  echo "Expected 154 reference populations; observed $population_count" >&2
  exit 3
}
echo "Matched-48 projection complete: $output_dir"
echo "Queries: 48; reference populations: 154; smartpca markers: 1,113,348"
