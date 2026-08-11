#!/usr/bin/env python3
"""Recalculate the 20 archived PCA recovery checks from public-safe inputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SUMMARY = ROOT / "data" / "summary" / "pca"
PROCESSED = ROOT / "data" / "processed" / "pca"
TITRATION = ROOT / "figures" / "supplementary_07_08_pca_titration" / "inputs" / "query_coordinates.tsv"
RECOVERY = ROOT / "workflows" / "10_pca" / "recovery"


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metric_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], float]:
    return {
        (row["metric"], row["space"], row["representation"]): float(row["value"])
        for row in rows
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raster-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    recomputed_metrics = args.output.parent / "pca_metrics_recomputed.csv"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflows" / "10_pca" / "verify_pca.py"),
            str(SUMMARY / "pca_ort_queries_matched48.evec"),
            str(SUMMARY / "pca_aadr_matched48.eval"),
            str(SUMMARY / "query_manifest_matched48.tsv"),
            str(SUMMARY / "reference_population_aggregates.tsv"),
            str(recomputed_metrics),
            "--reference",
            str(SUMMARY / "pca_metrics_matched48.csv"),
            "--absolute-tolerance",
            str(args.absolute_tolerance),
        ],
        check=True,
    )

    metrics = metric_lookup(read_csv(recomputed_metrics))
    build_summary = {
        row["metric"]: row["value"]
        for row in read_csv(PROCESSED / "projection_build_summary.tsv", delimiter="\t")
    }
    eigenvalues = [float(value) for value in (SUMMARY / "pca_aadr_matched48.eval").read_text().split()]
    populations = [
        line for line in (RECOVERY / "keep_pops_ref.txt").read_text().splitlines() if line
    ]
    query_rows = read_csv(SUMMARY / "query_manifest_matched48.tsv", delimiter="\t")
    titration_rows = read_csv(TITRATION, delimiter="\t")
    titration_ids = {row["iid"] for row in titration_rows}

    results: list[dict[str, str]] = []
    failures: list[str] = []

    def exact(check: str, manuscript: int, reproduced: int) -> None:
        passed = manuscript == reproduced
        results.append(
            {
                "check": check,
                "manuscript_value": str(manuscript),
                "reproduced_value": str(reproduced),
                "absolute_difference": str(abs(manuscript - reproduced)),
                "status": "PASS" if passed else "FAIL",
            }
        )
        if not passed:
            failures.append(check)

    def rounded(check: str, manuscript: float, reproduced: float, tolerance: float) -> None:
        difference = abs(manuscript - reproduced)
        passed = difference <= tolerance
        results.append(
            {
                "check": check,
                "manuscript_value": str(manuscript),
                "reproduced_value": repr(reproduced),
                "absolute_difference": repr(difference),
                "status": "PASS (rounding)" if passed else "FAIL",
            }
        )
        if not passed:
            failures.append(check)

    exact("reference_individuals", 985, len(eigenvalues))
    exact("reference_populations", 154, len(populations))
    exact("projected_query_points", 48, len(query_rows))
    exact("smartpca_markers", 1_113_348, int(build_summary["smartpca_markers_used"]))

    for index, manuscript in enumerate((6.45, 4.37, 1.53, 1.04), start=1):
        rounded(
            f"PC{index}_variance_pct",
            manuscript,
            metrics[("variance_explained_pct", "eigen", f"PC{index}")],
            0.005,
        )

    rounded(
        "imputed_spread_median_PC1_2",
        0.00047,
        metrics[("within_group_spread_median", "PC1-2", "imputed")],
        0.00001,
    )
    rounded(
        "imputed_spread_maximum_PC1_2",
        0.00102,
        metrics[("within_group_spread_maximum", "PC1-2", "imputed")],
        0.00001,
    )
    rounded(
        "pseudohaploid_spread_median_PC1_2",
        0.00217,
        metrics[("within_group_spread_median", "PC1-2", "pseudohaploid")],
        0.00001,
    )
    rounded(
        "pseudohaploid_spread_maximum_PC1_2",
        0.00336,
        metrics[("within_group_spread_maximum", "PC1-2", "pseudohaploid")],
        0.00001,
    )
    rounded(
        "imputed_vs_pseudohaploid_median_PC1_2",
        0.00309,
        metrics[("imputed_vs_pseudohaploid_median", "PC1-2", "matched24")],
        0.00001,
    )

    exact("titration_projected_points", 1694, len(titration_rows))
    exact(
        "titration_ORT15_points",
        842,
        sum(row["iid"].startswith("ORT15_") for row in titration_rows),
    )
    exact(
        "titration_ORT16_points",
        852,
        sum(row["iid"].startswith("ORT16_") for row in titration_rows),
    )
    missing_titration = sum(
        f"ORT15_nu_trim10__c0.2__r{replicate:02d}" not in titration_ids
        for replicate in range(1, 11)
    )
    exact("titration_missing_ORT15_nu_trim10_0.20x", 10, missing_titration)

    raster_rows = read_csv(PROCESSED / "raster_reference_checksums.tsv", delimiter="\t")
    raster_checks = (
        ("main_Figure_7_150dpi_raster", raster_rows[0]),
        ("Supplementary_Figure_S7_150dpi_raster", raster_rows[1]),
        ("Supplementary_Figure_S8_150dpi_raster", raster_rows[2]),
    )
    for check, expected in raster_checks:
        raster = args.raster_dir / expected["raster_filename"]
        observed = sha256(raster) if raster.is_file() else "missing"
        passed = observed == expected["sha256"]
        results.append(
            {
                "check": check,
                "manuscript_value": "identical",
                "reproduced_value": "identical" if passed else observed,
                "absolute_difference": "0 differing bytes" if passed else "checksum mismatch",
                "status": "PASS" if passed else "FAIL",
            }
        )
        if not passed:
            failures.append(check)

    archived = read_csv(PROCESSED / "validation_results.tsv", delimiter="\t")
    archived_checks = [row["check"] for row in archived]
    observed_checks = [row["check"] for row in results]
    if archived_checks != observed_checks:
        failures.append("validation-check inventory")

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "check",
                "manuscript_value",
                "reproduced_value",
                "absolute_difference",
                "status",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)

    if failures:
        raise SystemExit("PCA validation failed: " + ", ".join(failures))
    if len(results) != 20:
        raise SystemExit(f"Expected 20 checks; produced {len(results)}")
    print(f"PASS: {len(results)}/{len(results)} PCA recovery checks")
    print(f"PASS: metric-table comparison absolute tolerance {args.absolute_tolerance:g}")


if __name__ == "__main__":
    main()
