#!/usr/bin/env python3
"""Verify repository-source status for all seven main figures and Table 1."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "docs/main_output_source_manifest.tsv"
EXPECTED = {f"Figure {number}" for number in range(1, 8)} | {"Table 1"}
EXTERNAL_AUTHOR_ARTWORK = {"Figure 1"}


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
    repository_backed = 0
    for row in rows:
        if row["manuscript_output"] in EXTERNAL_AUTHOR_ARTWORK:
            for field in (
                "source_data_file",
                "script_used",
                "expected_output",
                "repository_path",
            ):
                if row[field].strip() != "NA":
                    raise AssertionError(
                        f"{row['manuscript_output']} must mark {field} as NA"
                    )
            if not row["validation_status"].strip():
                raise AssertionError(
                    f"missing validation status: {row['manuscript_output']}"
                )
            checked += 1
            continue

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
        repository_backed += 1

    print(
        "PASS  main-output source/status coverage: "
        f"{checked}/8 ({repository_backed} repository-backed; "
        "Figure 1 author-prepared externally)"
    )


if __name__ == "__main__":
    main()
