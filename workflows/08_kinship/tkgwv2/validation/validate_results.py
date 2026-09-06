#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Validate all 16 TKGWV2 manuscript comparisons exactly."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
from pathlib import Path


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ROOT = PACKAGE.parents[2]
EXPECTED_PATH = ROOT / "data" / "summary" / "kinship" / "tkgwv2_pair_results.tsv"
RECOVERED_VALIDATION_PATH = HERE / "reproduced_results.tsv"
CHECKED_REPORT_PATH = HERE / "validation_results.tsv"
RUNS = ("fu_fu", "nu_nu", "nu_fu")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_relationship(value: str) -> str:
    normalized = value.strip().lower().replace(" ", "_")
    return {
        "1st_degree": "first_degree",
        "2nd_degree": "second_degree",
    }.get(normalized, normalized)


def expected_run_pairs() -> set[tuple[str, str, str]]:
    pairs: set[tuple[str, str, str]] = set()
    for run in RUNS:
        path = PACKAGE / "dyads" / f"{run}.tsv"
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 2:
                    raise SystemExit(f"invalid dyad at {path.name}:{line_number}")
                pairs.add((run, fields[0], fields[1]))
    if len(pairs) != 16:
        raise SystemExit(f"dyad files must contain 16 unique run/pair rows, found {len(pairs)}")
    return pairs


def load_tkgwv2_outputs(results_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for run in RUNS:
        path = results_dir / run / "TKGWV2_Results.txt"
        if not path.is_file():
            raise SystemExit(f"missing TKGWV2 result: {run}/TKGWV2_Results.txt")
        for raw in read_tsv(path):
            rows.append(
                {
                    "run": run,
                    "sample1": raw["Sample1"],
                    "sample2": raw["Sample2"],
                    "used_snps": raw["Used_SNPs"],
                    "HRC": raw["HRC"],
                    "counts0": raw["counts0"],
                    "counts4": raw["counts4"],
                    "classification": normalize_relationship(raw["Relationship"]),
                    "fresh_matches_recovered_fields": "NOT_CHECKED",
                }
            )
    return rows


def build_report(
    expected_rows: list[dict[str, str]], observed_rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    expected = {(r["sample1"], r["sample2"]): r for r in expected_rows}
    observed = {(r["sample1"], r["sample2"]): r for r in observed_rows}
    if len(expected_rows) != 16 or len(expected) != 16:
        raise SystemExit("expected-results table must contain 16 unique pairs")
    if len(observed_rows) != 16 or len(observed) != 16:
        raise SystemExit("reproduced results must contain 16 unique pairs")
    if set(expected) != set(observed):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise SystemExit(f"pair-set mismatch; missing={missing}, extra={extra}")
    observed_run_pairs = {
        (row["run"], row["sample1"], row["sample2"]) for row in observed_rows
    }
    if observed_run_pairs != expected_run_pairs():
        raise SystemExit("run/dyad membership does not match the recovered definitions")

    report: list[dict[str, str]] = []
    for row in observed_rows:
        pair = (row["sample1"], row["sample2"])
        target = expected[pair]
        used_diff = abs(int(row["used_snps"]) - int(target["used_snps"]))
        hrc_diff = abs(Decimal(row["HRC"]) - Decimal(target["HRC"]))
        counts0_diff = abs(int(row["counts0"]) - int(target["counts0"]))
        counts4_diff = abs(int(row["counts4"]) - int(target["counts4"]))
        classification_match = (
            normalize_relationship(row["classification"])
            == normalize_relationship(target["classification"])
        )
        scientific_match = (
            used_diff == 0
            and hrc_diff == 0
            and counts0_diff == 0
            and counts4_diff == 0
            and classification_match
        )
        source_match = row.get("fresh_matches_recovered_fields", "NOT_CHECKED")
        provenance_match = source_match in {"PASS", "NOT_CHECKED"}
        report.append(
            {
                "run": row["run"],
                "sample1": pair[0],
                "sample2": pair[1],
                "expected_used_snps": target["used_snps"],
                "reproduced_used_snps": row["used_snps"],
                "used_snps_abs_diff": str(used_diff),
                "expected_hrc": target["HRC"],
                "reproduced_hrc": row["HRC"],
                "hrc_abs_diff": f"{hrc_diff:.4f}",
                "expected_counts0": target["counts0"],
                "reproduced_counts0": row["counts0"],
                "counts0_abs_diff": str(counts0_diff),
                "expected_counts4": target["counts4"],
                "reproduced_counts4": row["counts4"],
                "counts4_abs_diff": str(counts4_diff),
                "classification_match": "PASS" if classification_match else "FAIL",
                "fresh_matches_recovered_fields": source_match,
                "status": "PASS" if scientific_match and provenance_match else "FAIL",
            }
        )

    focal = observed[("ORT15_fu_bamrefine5", "ORT16_fu_bamrefine5")]
    if focal["HRC"] != "0.2356" or focal["used_snps"] != "4517214":
        raise SystemExit("focal HRC or used-SNP count does not match the manuscript")
    return report


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--results-dir",
        type=Path,
        help="directory containing fu_fu, nu_nu and nu_fu TKGWV2 outputs",
    )
    source.add_argument(
        "--reproduced-results",
        type=Path,
        default=RECOVERED_VALIDATION_PATH,
        help="normalized reproduced-results TSV (default: checked-in recovery output)",
    )
    parser.add_argument("--output", type=Path, help="write the comparison TSV here")
    args = parser.parse_args()

    observed = (
        load_tkgwv2_outputs(args.results_dir)
        if args.results_dir
        else read_tsv(args.reproduced_results)
    )
    report = build_report(read_tsv(EXPECTED_PATH), observed)

    if args.output:
        write_report(args.output, report)
    elif not args.results_dir and args.reproduced_results == RECOVERED_VALIDATION_PATH:
        if report != read_tsv(CHECKED_REPORT_PATH):
            raise SystemExit("checked validation_results.tsv is stale")

    failures = [row for row in report if row["status"] != "PASS"]
    print(f"TKGWV2 validation: {len(report) - len(failures)}/16 PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
