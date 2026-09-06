#!/usr/bin/env python3
"""Validate the canonical collapsed-only mitochondrial BAM manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


EXPECTED_FIELDS = ["sample", "read_set_state", "input_scope", "bam"]
EXPECTED_STATE = "full_udg_collapsed_only"
EXPECTED_SCOPE = "combined_libraries_a_b_c"
NESTED_FILENAME = re.compile(r"[.](?:2|3)(?:[._-]|$)", re.IGNORECASE)
KNOWN_CHECKSUMS = Path(__file__).resolve().parent / "corrected_input_checksums.tsv"


def known_checksums(path: Path = KNOWN_CHECKSUMS) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["sample", "read_set_state", "input_scope", "sha256"]
        if reader.fieldnames != expected:
            raise ValueError(f"unexpected checksum-table header: {reader.fieldnames}")
        rows = list(reader)
    result = {row["sample"]: row["sha256"] for row in rows}
    if set(result) != {"ORT15", "ORT16"} or any(
        row["read_set_state"] != EXPECTED_STATE
        or row["input_scope"] != EXPECTED_SCOPE
        or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"])
        for row in rows
    ):
        raise ValueError("checksum table does not define the verified canonical inputs")
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest(path: Path, require_files: bool = False) -> list[dict[str, str]]:
    expected_sha = known_checksums()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != EXPECTED_FIELDS:
            raise ValueError(
                f"expected manifest header {EXPECTED_FIELDS}, found {reader.fieldnames}"
            )
        rows = []
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise ValueError(f"line {line_number}: too many tab-separated fields")
            row = {key: (value or "").strip() for key, value in raw_row.items()}
            if any(row.values()):
                rows.append(row)

    if not rows:
        raise ValueError("manifest has no data rows")

    samples: set[str] = set()
    bams: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        sample = row["sample"]
        bam = row["bam"]
        if not sample or not bam:
            raise ValueError(f"line {line_number}: sample and BAM are required")
        if sample not in expected_sha:
            raise ValueError(f"line {line_number}: unsupported sample: {sample}")
        if row["read_set_state"] != EXPECTED_STATE:
            raise ValueError(
                f"line {line_number}: only {EXPECTED_STATE!r} is permitted; "
                "nested read-set supersets must not be merged"
            )
        if row["input_scope"] != EXPECTED_SCOPE:
            raise ValueError(
                f"line {line_number}: input_scope must be {EXPECTED_SCOPE!r}"
            )
        if NESTED_FILENAME.search(Path(bam).name):
            raise ValueError(
                f"line {line_number}: .2/.3 nested-state BAMs are prohibited: {bam}"
            )
        if sample in samples:
            raise ValueError(f"line {line_number}: duplicate sample row: {sample}")
        if bam in bams:
            raise ValueError(f"line {line_number}: BAM reused by more than one sample: {bam}")
        if require_files:
            bam_path = Path(bam)
            if not bam_path.is_file():
                raise ValueError(f"line {line_number}: BAM does not exist: {bam}")
            if sha256(bam_path) != expected_sha[sample]:
                raise ValueError(
                    f"line {line_number}: BAM SHA-256 is not the verified {sample} "
                    "collapsed-only input"
                )
        samples.add(sample)
        bams.add(bam)

    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--sample")
    parser.add_argument("--require-files", action="store_true")
    args = parser.parse_args()

    rows = validate_manifest(args.manifest, require_files=args.require_files)
    if args.sample:
        matches = [row for row in rows if row["sample"] == args.sample]
        if len(matches) != 1:
            raise ValueError(f"expected one manifest row for {args.sample}, found {len(matches)}")
        print(matches[0]["bam"])
    else:
        print(f"PASS  canonical mtDNA BAM manifest: {len(rows)} unique sample rows")


if __name__ == "__main__":
    main()
