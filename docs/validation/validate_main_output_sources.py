#!/usr/bin/env python3
"""Verify that all seven main figures and Table 1 have identifiable sources."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "docs/main_output_source_manifest.tsv"
EXPECTED = {f"Figure {number}" for number in range(1, 8)} | {"Table 1"}


def split_paths(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def main() -> None:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    observed = {row["manuscript_output"] for row in rows}
    if len(rows) != 8 or observed != EXPECTED:
        raise AssertionError(
            f"expected Figure 1-7 and Table 1, found {sorted(observed)}"
        )

    checked = 0
    for row in rows:
        for field in ("source_data_file", "script_used"):
            paths = split_paths(row[field])
            if not paths:
                raise AssertionError(
                    f"{row['manuscript_output']} has no {field.replace('_', ' ')}"
                )
            for relative_path in paths:
                path = ROOT / relative_path
                if not path.is_file():
                    raise AssertionError(f"missing mapped file: {relative_path}")
        repository_path = ROOT / row["repository_path"]
        if not repository_path.exists():
            raise AssertionError(
                f"missing repository path: {row['repository_path']}"
            )
        if not row["validation_status"].strip():
            raise AssertionError(
                f"missing validation status: {row['manuscript_output']}"
            )
        checked += 1

    print(f"PASS  main-output source coverage: {checked}/8")


if __name__ == "__main__":
    main()
