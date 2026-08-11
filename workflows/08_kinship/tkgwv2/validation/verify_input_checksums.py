#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify the 32 PED/MAP files and EAS frequency file for exact reruns."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parent.parent
MANIFEST = PACKAGE / "input_checksums.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ped-map-dir", required=True, type=Path)
    parser.add_argument("--frequency-file", required=True, type=Path)
    args = parser.parse_args()

    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 33:
        raise SystemExit(f"input manifest must contain 33 rows, found {len(rows)}")

    failures: list[str] = []
    for row in rows:
        path = (
            args.ped_map_dir / row["file"]
            if row["class"] == "processed_input"
            else args.frequency_file
        )
        if not path.is_file():
            failures.append(f"{row['file']}: missing")
            continue
        if path.stat().st_size != int(row["bytes"]):
            failures.append(f"{row['file']}: byte-size mismatch")
            continue
        if sha256(path) != row["sha256"]:
            failures.append(f"{row['file']}: SHA-256 mismatch")

    if failures:
        for failure in failures:
            print(f"FAIL\t{failure}")
        return 1
    print("TKGWV2 input verification: 33/33 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
