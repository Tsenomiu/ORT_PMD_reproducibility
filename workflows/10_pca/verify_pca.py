#!/usr/bin/env python3
"""Verify released ORT-only PCA coordinates without reference-individual rows."""
from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


IDENTITY_FIELDS = ("metric", "space", "individual", "library", "representation", "n_points")


def max_pairwise(points: list[list[float]]) -> float:
    return max((math.dist(points[i], points[j]) for i in range(len(points)) for j in range(i + 1, len(points))), default=float("nan"))


def compare_metric_tables(reference: Path, observed: Path, abs_tol: float) -> None:
    with reference.open(newline="") as handle:
        expected_rows = list(csv.DictReader(handle))
    with observed.open(newline="") as handle:
        observed_rows = list(csv.DictReader(handle))
    if len(expected_rows) != len(observed_rows):
        raise ValueError(
            f"PCA metric row count differs: expected {len(expected_rows)}, observed {len(observed_rows)}"
        )
    for row_number, (expected, observed) in enumerate(zip(expected_rows, observed_rows), start=2):
        for field in IDENTITY_FIELDS:
            if expected[field] != observed[field]:
                raise ValueError(
                    f"PCA metric identity differs at CSV row {row_number}, field {field}: "
                    f"expected {expected[field]!r}, observed {observed[field]!r}"
                )
        expected_value = float(expected["value"])
        observed_value = float(observed["value"])
        if not math.isclose(expected_value, observed_value, rel_tol=0.0, abs_tol=abs_tol):
            raise ValueError(
                f"PCA metric value differs at CSV row {row_number}: "
                f"expected {expected_value:.17g}, observed {observed_value:.17g}, "
                f"absolute tolerance {abs_tol:g}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("evec", type=Path); parser.add_argument("eval", type=Path)
    parser.add_argument("manifest", type=Path); parser.add_argument("reference_aggregates", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference", type=Path, help="Archived metric table to compare with the reconstructed table")
    parser.add_argument("--absolute-tolerance", type=float, default=1e-12,
                        help="Maximum absolute floating-point difference allowed with --reference (default: 1e-12)")
    args = parser.parse_args()
    values = [float(x) for x in args.eval.read_text().split()]; trace = sum(x for x in values if x > 0)
    if len(values) != 985 or abs(trace - 984.0) > 0.001: raise ValueError("Expected 985 eigenvalues and full trace 984")
    with args.manifest.open(newline="") as handle:
        manifest = {r["prefixed_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with args.reference_aggregates.open(newline="") as handle:
        aggregate_rows = {r["population"]: r for r in csv.DictReader(handle, delimiter="\t")}
    yakut_aggregate = aggregate_rows["Yakut.DG"]
    groups = defaultdict(list); queries = {}; imp = phap = 0
    with args.evec.open() as handle:
        next(handle)
        for raw in handle:
            fields = raw.split(); iid = fields[0]; pcs = list(map(float, fields[1:11]))
            if iid in manifest:
                item = manifest[iid]; groups[(item["individual"], item["library"], item["representation"])].append(pcs)
                queries[iid] = pcs
                imp += item["representation"] == "imputed"; phap += item["representation"] == "pseudohaploid"
            else:
                raise ValueError(f"Non-query individual row must not be in the public evec: {iid}")
    if (imp, phap) != (24, 24): raise ValueError(f"Unexpected query counts: {(imp, phap)}")
    rows = []
    for pc_count, space in ((2, "PC1-2"), (4, "PC1-4")):
        spreads_by_representation = defaultdict(list)
        for key, points in sorted(groups.items()):
            spread = max_pairwise([p[:pc_count] for p in points])
            rows.append({"metric": "within_group_spread", "space": space, "individual": key[0], "library": key[1], "representation": key[2],
                         "value": spread, "n_points": len(points)})
            spreads_by_representation[key[2]].append(spread)
        for representation, spreads in sorted(spreads_by_representation.items()):
            rows.append({"metric": "within_group_spread_median", "space": space, "individual": "", "library": "", "representation": representation,
                         "value": statistics.median(spreads), "n_points": len(spreads)})
            rows.append({"metric": "within_group_spread_maximum", "space": space, "individual": "", "library": "", "representation": representation,
                         "value": max(spreads), "n_points": len(spreads)})
        matched = []
        for iid, item in manifest.items():
            if item["representation"] != "imputed": continue
            other = "PHAP__" + iid.split("__", 1)[1]
            matched.append(math.dist(queries[iid][:pc_count], queries[other][:pc_count]))
        rows.append({"metric": "imputed_vs_pseudohaploid_median", "space": space, "individual": "", "library": "", "representation": "matched24",
                     "value": statistics.median(matched), "n_points": len(matched)})
        all_queries = list(queries.values())
        qcentroid = [sum(p[i] for p in all_queries) / len(all_queries) for i in range(pc_count)]
        ycentroid = [float(yakut_aggregate[f"PC{i+1}_mean"]) for i in range(pc_count)]
        ort_yakut = math.dist(qcentroid, ycentroid)
        yakut_scatter = float(yakut_aggregate[f"mean_distance_to_centroid_PC1_{pc_count}"])
        yakut_n = int(yakut_aggregate["n"])
        rows.append({"metric": "ORT_to_Yakut_centroid_distance", "space": space, "individual": "", "library": "", "representation": "all48",
                     "value": ort_yakut, "n_points": len(all_queries)})
        rows.append({"metric": "within_Yakut_mean_to_centroid", "space": space, "individual": "", "library": "", "representation": "reference",
                     "value": yakut_scatter, "n_points": yakut_n})
    for index in range(4):
        rows.append({"metric": "variance_explained_pct", "space": "eigen", "individual": "", "library": "", "representation": f"PC{index+1}",
                     "value": 100 * values[index] / trace, "n_points": 1})
    headline = {(r["metric"], r["representation"]): float(r["value"]) for r in rows if r["space"] == "PC1-2"}
    checks = {
        ("within_group_spread_median", "imputed"): 0.00047,
        ("within_group_spread_maximum", "imputed"): 0.00102,
        ("within_group_spread_median", "pseudohaploid"): 0.00217,
        ("within_group_spread_maximum", "pseudohaploid"): 0.00336,
        ("imputed_vs_pseudohaploid_median", "matched24"): 0.00309,
        ("ORT_to_Yakut_centroid_distance", "all48"): 0.00532,
        ("within_Yakut_mean_to_centroid", "reference"): 0.00181,
    }
    for key, expected in checks.items():
        if abs(headline[key] - expected) > 0.00001:
            raise ValueError(f"Headline PCA metric {key}={headline[key]:.8f}, expected about {expected:.5f}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    if args.reference is not None:
        compare_metric_tables(args.reference, args.output, args.absolute_tolerance)


if __name__ == "__main__":
    main()
