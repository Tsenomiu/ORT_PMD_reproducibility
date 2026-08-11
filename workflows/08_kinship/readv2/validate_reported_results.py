#!/usr/bin/env python3
"""Validate the curated READv2 derived-data package.

Generated publication validator; this is not part of upstream READv2.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[3]
DATA = REPOSITORY / "data" / "processed" / "readv2"


EXPECTED = {
    "primary_chr1_22_X_Y": {
        "p0_file": "primary_p0_anonymized.tsv",
        "pair_count": 210,
        "baseline": 0.24671998039541876,
        "raw_p0": 0.18973589010001654,
        "p0_norm": 0.7690333381022743,
        "kc": 0.23096666189772574,
        "overlap": 169172,
        "degree": "First Degree",
        "subtype": "Parent-offspring",
        "rounded_p0": "0.769",
        "rounded_kc": "0.231",
    },
    "autosome_1240k": {
        "p0_file": "autosome_1240k_p0_anonymized.tsv",
        "pair_count": 230,
        "baseline": 0.24853168423505886,
        "raw_p0": 0.19097369975420925,
        "p0_norm": 0.7684078605188551,
        "kc": 0.23159213948114488,
        "overlap": 157451,
        "degree": "First Degree",
        "subtype": "Parent-offspring",
        "rounded_p0": "0.768",
        "rounded_kc": "0.232",
    },
    "transversion_only": {
        "p0_file": "transversion_only_p0_anonymized.tsv",
        "pair_count": 230,
        "baseline": 0.24893068550585545,
        "raw_p0": 0.19167202572347267,
        "p0_norm": 0.7699815124598774,
        "kc": 0.23001848754012255,
        "overlap": 31100,
        "degree": "First Degree",
        "subtype": "N/A",
        "rounded_p0": "0.770",
        "rounded_kc": "0.230",
    },
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


checks: list[dict[str, str]] = []


def record(analysis: str, check: str, observed: object, expected: object, passed: bool) -> None:
    checks.append(
        {
            "analysis": analysis,
            "check": check,
            "recomputed_value": str(observed),
            "reported_value": str(expected),
            "status": "PASS" if passed else "FAIL",
        }
    )


def close(observed: float, expected: float) -> bool:
    return math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-15)


def readv2_call(p0_norm: float, window_fraction: float, effective_snps: float) -> tuple[str, str]:
    degree = "First Degree" if 0.625 <= p0_norm <= 0.8125 else "not first degree"
    subtype = "N/A"
    if degree == "First Degree" and effective_snps >= 10000:
        if window_fraction >= 0.6:
            subtype = "NotApplicable"
        elif window_fraction >= 0.35:
            subtype = "Siblings"
        elif window_fraction <= 0.3:
            subtype = "Parent-offspring"
    return degree, subtype


def validate_samples() -> None:
    rows = read_tsv(DATA / "sample_manifest.tsv")
    primary = [row for row in rows if row["analysis"] == "primary_chr1_22_X_Y"]
    later = [row for row in rows if row["analysis"] == "autosome_and_transversion"]
    removed = sorted(
        row["sample_id"] for row in primary if row["status"] == "removed_by_mind_0.8"
    )
    retained = [row for row in primary if row["status"] == "retained_after_mind_0.8"]
    expected_removed = ["ORT1", "ORT17", "ORT19", "ORT20", "ORT22", "ORT25"]
    record("sample_manifest", "primary_input_people", len(primary), 27, len(primary) == 27)
    record("sample_manifest", "primary_retained_people", len(retained), 21, len(retained) == 21)
    record(
        "sample_manifest",
        "primary_mind_exclusions",
        ",".join(removed),
        ",".join(expected_removed),
        removed == expected_removed,
    )
    record("sample_manifest", "later_input_people", len(later), 22, len(later) == 22)


def validate_results() -> None:
    focal = {row["analysis"]: row for row in read_tsv(DATA / "focal_results.tsv")}
    summary = {row["analysis"]: row for row in read_tsv(DATA / "analysis_summary.tsv")}

    for analysis, expected in EXPECTED.items():
        pairwise = read_tsv(DATA / str(expected["p0_file"]))
        ranks = [int(row["rank"]) for row in pairwise]
        raw_values = [float(row["nonnormalized_p0"]) for row in pairwise]
        normalized_values = [float(row["normalized_p0"]) for row in pairwise]
        baseline = statistics.median(raw_values)

        record(
            analysis,
            "pair_count",
            len(pairwise),
            expected["pair_count"],
            len(pairwise) == expected["pair_count"],
        )
        record(
            analysis,
            "rank_sequence",
            f"1-{len(ranks)}",
            f"1-{expected['pair_count']}",
            ranks == list(range(1, len(ranks) + 1)),
        )
        record(
            analysis,
            "all_pair_normalized_values",
            "raw_P0 / cohort_median",
            "raw_P0 / cohort_median",
            all(close(norm, raw / baseline) for raw, norm in zip(raw_values, normalized_values)),
        )
        record(
            analysis,
            "median_baseline",
            format(baseline, ".17g"),
            format(float(expected["baseline"]), ".17g"),
            close(baseline, float(expected["baseline"])),
        )

        row = focal[analysis]
        raw_p0 = float(row["Nonnormalized_P0"])
        p0_norm = raw_p0 / baseline
        kc = 1.0 - p0_norm
        degree, subtype = readv2_call(
            p0_norm,
            float(row["Perc_Win_1stdeg_P0"]),
            float(row["NSNPsXNorm"]),
        )

        record(
            analysis,
            "focal_raw_p0",
            format(raw_p0, ".17g"),
            format(float(expected["raw_p0"]), ".17g"),
            close(raw_p0, float(expected["raw_p0"])),
        )
        record(
            analysis,
            "normalized_p0",
            format(p0_norm, ".17g"),
            row["P0_mean"],
            close(p0_norm, float(row["P0_mean"]))
            and close(p0_norm, float(expected["p0_norm"])),
        )
        record(
            analysis,
            "kinship_coefficient",
            format(kc, ".17g"),
            row["KinshipCoefficient"],
            close(kc, float(row["KinshipCoefficient"]))
            and close(kc, float(expected["kc"])),
        )
        record(
            analysis,
            "overlap_snps",
            int(row["OverlapNSNPs"]),
            expected["overlap"],
            int(row["OverlapNSNPs"]) == expected["overlap"],
        )
        record(analysis, "degree", degree, expected["degree"], degree == expected["degree"])
        record(analysis, "subtype", subtype, expected["subtype"], subtype == expected["subtype"])
        record(
            analysis,
            "rounded_p0_3dp",
            f"{p0_norm:.3f}",
            expected["rounded_p0"],
            f"{p0_norm:.3f}" == expected["rounded_p0"],
        )
        record(
            analysis,
            "rounded_kc_3dp",
            f"{kc:.3f}",
            expected["rounded_kc"],
            f"{kc:.3f}" == expected["rounded_kc"],
        )

        summary_row = summary[analysis]
        summary_ok = (
            close(float(summary_row["median_raw_p0"]), baseline)
            and close(float(summary_row["normalized_p0"]), p0_norm)
            and close(float(summary_row["kinship_coefficient"]), kc)
            and int(summary_row["overlap_snps"]) == int(expected["overlap"])
            and summary_row["degree"] == degree
            and summary_row["subtype"] == subtype
        )
        record(analysis, "analysis_summary", "consistent", "consistent", summary_ok)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DATA / "validation_results.tsv",
        help="Machine-readable validation output (default: tracked validation_results.tsv)",
    )
    args = parser.parse_args()

    validate_samples()
    validate_results()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["analysis", "check", "recomputed_value", "reported_value", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(checks)

    failures = [row for row in checks if row["status"] != "PASS"]
    if failures:
        raise SystemExit(f"FAIL: {len(failures)} of {len(checks)} READv2 checks failed")
    print(f"PASS: {len(checks)}/{len(checks)} READv2 checks; wrote {args.output}")


if __name__ == "__main__":
    main()
