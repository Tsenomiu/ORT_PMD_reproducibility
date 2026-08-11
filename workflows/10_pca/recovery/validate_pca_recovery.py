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

import numpy as np
from PIL import Image, ImageFilter


ROOT = Path(__file__).resolve().parents[3]
SUMMARY = ROOT / "data" / "summary" / "pca"
PROCESSED = ROOT / "data" / "processed" / "pca"
TITRATION = ROOT / "figures" / "supplementary_07_08_pca_titration" / "inputs" / "query_coordinates.tsv"
RECOVERY = ROOT / "workflows" / "10_pca" / "recovery"
CANONICAL_RASTERS = RECOVERY / "canonical_rasters"

# These thresholds tolerate renderer-version anti-aliasing differences while
# rejecting missing, displaced or visibly altered plot content. The canonical
# strict mode remains the default and continues to require exact SHA-256 values.
PORTABLE_MAX_NMAE = 0.020
PORTABLE_MAX_DHASH_FRACTION = 0.060
PORTABLE_MIN_PIXEL_CORRELATION = 0.970
PORTABLE_MIN_REFERENCE_INK_COVERAGE = 0.990
PORTABLE_MIN_OBSERVED_INK_COVERAGE = 0.998
PORTABLE_INK_THRESHOLD = 245
PORTABLE_DILATION_SIZE = 5


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def difference_hash(image: Image.Image, size: int = 32) -> np.ndarray:
    gray = image.convert("L").resize((size + 1, size), Image.Resampling.LANCZOS)
    pixels = np.asarray(gray, dtype=np.int16)
    return (pixels[:, 1:] > pixels[:, :-1]).reshape(-1)


def ink_mask(image: Image.Image) -> np.ndarray:
    gray = np.asarray(image.convert("L"), dtype=np.uint8)
    return gray < PORTABLE_INK_THRESHOLD


def dilate(mask: np.ndarray) -> np.ndarray:
    image = Image.fromarray(mask.astype(np.uint8) * 255)
    expanded = image.filter(ImageFilter.MaxFilter(PORTABLE_DILATION_SIZE))
    return np.asarray(expanded, dtype=np.uint8) > 0


def coverage(source: np.ndarray, target_dilated: np.ndarray) -> float:
    denominator = int(source.sum())
    if denominator == 0:
        return 1.0 if int(target_dilated.sum()) == 0 else 0.0
    return float(np.logical_and(source, target_dilated).sum() / denominator)


def portable_raster_metrics(reference_path: Path, observed_path: Path) -> dict[str, object]:
    with Image.open(reference_path) as reference_file:
        reference_format = reference_file.format
        reference_mode = reference_file.mode
        reference = reference_file.convert("RGB")
    with Image.open(observed_path) as observed_file:
        observed_format = observed_file.format
        observed_mode = observed_file.mode
        observed = observed_file.convert("RGB")

    metrics: dict[str, object] = {
        "reference_format": reference_format,
        "observed_format": observed_format,
        "reference_mode": reference_mode,
        "observed_mode": observed_mode,
        "reference_size": reference.size,
        "observed_size": observed.size,
    }
    if reference.size != observed.size:
        return metrics

    reference_pixels = np.asarray(reference, dtype=np.float32)
    observed_pixels = np.asarray(observed, dtype=np.float32)
    metrics["nmae"] = float(np.mean(np.abs(reference_pixels - observed_pixels)) / 255.0)
    correlation_size = (
        max(1, reference.width // 4),
        max(1, reference.height // 4),
    )
    reference_correlation = np.asarray(
        reference.convert("L").resize(correlation_size, Image.Resampling.LANCZOS),
        dtype=np.float32,
    )
    observed_correlation = np.asarray(
        observed.convert("L").resize(correlation_size, Image.Resampling.LANCZOS),
        dtype=np.float32,
    )
    reference_centered = reference_correlation.reshape(-1) - float(reference_correlation.mean())
    observed_centered = observed_correlation.reshape(-1) - float(observed_correlation.mean())
    correlation_denominator = float(
        np.linalg.norm(reference_centered) * np.linalg.norm(observed_centered)
    )
    metrics["pixel_correlation"] = (
        float(np.dot(reference_centered, observed_centered) / correlation_denominator)
        if correlation_denominator
        else 1.0 if np.array_equal(reference_pixels, observed_pixels) else 0.0
    )

    reference_hash = difference_hash(reference)
    observed_hash = difference_hash(observed)
    metrics["dhash_fraction"] = float(np.mean(reference_hash != observed_hash))

    reference_ink = ink_mask(reference)
    observed_ink = ink_mask(observed)
    metrics["reference_ink_coverage"] = coverage(reference_ink, dilate(observed_ink))
    metrics["observed_ink_coverage"] = coverage(observed_ink, dilate(reference_ink))
    return metrics


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
    parser.add_argument(
        "--mode",
        choices=("canonical-strict", "ci-portable"),
        default="canonical-strict",
        help="canonical-strict requires exact raster SHA-256; ci-portable uses documented visual-equivalence thresholds",
    )
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
        if args.mode == "canonical-strict":
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
        else:
            reference = CANONICAL_RASTERS / expected["raster_filename"]
            reference_valid = reference.is_file() and sha256(reference) == expected["sha256"]
            try:
                metrics = portable_raster_metrics(reference, raster)
            except (FileNotFoundError, OSError, ValueError) as error:
                metrics = {"error": str(error)}

            same_size = metrics.get("reference_size") == metrics.get("observed_size")
            png_format = metrics.get("reference_format") == metrics.get("observed_format") == "PNG"
            rgb_mode = metrics.get("reference_mode") == metrics.get("observed_mode") == "RGB"
            nmae = float(metrics.get("nmae", math.inf))
            pixel_correlation = float(metrics.get("pixel_correlation", -math.inf))
            dhash_fraction = float(metrics.get("dhash_fraction", math.inf))
            reference_coverage = float(metrics.get("reference_ink_coverage", -math.inf))
            observed_coverage = float(metrics.get("observed_ink_coverage", -math.inf))
            passed = (
                reference_valid
                and same_size
                and png_format
                and rgb_mode
                and nmae <= PORTABLE_MAX_NMAE
                and pixel_correlation >= PORTABLE_MIN_PIXEL_CORRELATION
                and dhash_fraction <= PORTABLE_MAX_DHASH_FRACTION
                and reference_coverage >= PORTABLE_MIN_REFERENCE_INK_COVERAGE
                and observed_coverage >= PORTABLE_MIN_OBSERVED_INK_COVERAGE
            )
            expected_size = metrics.get("reference_size", "missing")
            observed_size = metrics.get("observed_size", "missing")
            metric_summary = (
                f"{metrics.get('observed_format', 'missing')}/"
                f"{metrics.get('observed_mode', 'missing')} {observed_size}; "
                f"NMAE={nmae:.6f}; correlation={pixel_correlation:.6f}; "
                f"dHash={dhash_fraction:.6f}; "
                f"reference ink={reference_coverage:.6f}; observed ink={observed_coverage:.6f}"
            )
            print(f"CI-portable raster metrics [{check}]: {metric_summary}")
            results.append(
                {
                    "check": check,
                    "manuscript_value": (
                        f"canonical PNG {expected_size}; NMAE<={PORTABLE_MAX_NMAE:.3f}; "
                        f"correlation>={PORTABLE_MIN_PIXEL_CORRELATION:.3f}; "
                        f"dHash<={PORTABLE_MAX_DHASH_FRACTION:.3f}; "
                        f"reference ink>={PORTABLE_MIN_REFERENCE_INK_COVERAGE:.3f}; "
                        f"observed ink>={PORTABLE_MIN_OBSERVED_INK_COVERAGE:.3f}"
                    ),
                    "reproduced_value": (
                        metric_summary
                    ),
                    "absolute_difference": "platform-tolerant visual comparison",
                    "status": "PASS (ci-portable)" if passed else "FAIL",
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
    print(f"PASS: {len(results)}/{len(results)} PCA recovery checks ({args.mode})")
    print(f"PASS: metric-table comparison absolute tolerance {args.absolute_tolerance:g}")
    if args.mode == "ci-portable":
        print("PASS: three generated PNG rasters matched format, dimensions and portable visual-equivalence thresholds")


if __name__ == "__main__":
    main()
