#!/usr/bin/env python3
"""Write aggregate target-site read support without publishing read names."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
OUTPUT_FIELDS = [
    "sample",
    "chrom",
    "position",
    "ref",
    "alt",
    "mapq_min",
    "baseq_min",
    "qualifying_ref_reads",
    "qualifying_alt_reads",
    "qualifying_other_reads",
    "distinct_alt_read_names",
    "alt_forward",
    "alt_reverse",
    "alt_fraction_ref_alt",
    "alt_mapq_min",
    "alt_mapq_max",
    "alt_baseq_min",
    "alt_baseq_max",
]


def base_at_ref_pos(
    sequence: str, quality: str, start: int, cigar: str, target: int
) -> tuple[str | None, int | None]:
    reference_position = start
    query_position = 0
    for length_text, operation in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if operation in "M=X":
            if reference_position <= target < reference_position + length:
                offset = query_position + target - reference_position
                return sequence[offset].upper(), ord(quality[offset]) - 33
            reference_position += length
            query_position += length
        elif operation in "IS":
            query_position += length
        elif operation in "DN":
            if reference_position <= target < reference_position + length:
                return "*", None
            reference_position += length
        elif operation in "HP":
            continue
        else:  # pragma: no cover - regular expression limits operations
            raise ValueError(f"unsupported CIGAR operation: {operation}")
    return None, None


def summarize_site(
    samtools: str,
    bam: Path,
    sample: str,
    chrom: str,
    position: int,
    ref: str,
    alt: str,
    mapq_min: int,
    baseq_min: int,
) -> dict[str, str | int]:
    command = [samtools, "view", "-q", str(mapq_min), str(bam), f"{chrom}:{position}-{position}"]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)

    passing: list[dict[str, str | int]] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        # Match the default mpileup exclusions: unmapped, secondary, QC-fail and duplicate.
        if flag & (0x4 | 0x100 | 0x200 | 0x400):
            continue
        base, base_quality = base_at_ref_pos(
            fields[9], fields[10], int(fields[3]), fields[5], position
        )
        if base in (None, "*") or base_quality is None or base_quality < baseq_min:
            continue
        passing.append(
            {
                "qname": fields[0],
                "strand": "reverse" if flag & 0x10 else "forward",
                "mapq": int(fields[4]),
                "base": base,
                "baseq": base_quality,
            }
        )

    ref_reads = [row for row in passing if row["base"] == ref]
    alt_reads = [row for row in passing if row["base"] == alt]
    other_reads = [row for row in passing if row["base"] not in (ref, alt)]
    denominator = len(ref_reads) + len(alt_reads)
    alt_mapq = [int(row["mapq"]) for row in alt_reads]
    alt_baseq = [int(row["baseq"]) for row in alt_reads]
    return {
        "sample": sample,
        "chrom": chrom,
        "position": position,
        "ref": ref,
        "alt": alt,
        "mapq_min": mapq_min,
        "baseq_min": baseq_min,
        "qualifying_ref_reads": len(ref_reads),
        "qualifying_alt_reads": len(alt_reads),
        "qualifying_other_reads": len(other_reads),
        "distinct_alt_read_names": len({str(row["qname"]) for row in alt_reads}),
        "alt_forward": sum(row["strand"] == "forward" for row in alt_reads),
        "alt_reverse": sum(row["strand"] == "reverse" for row in alt_reads),
        "alt_fraction_ref_alt": len(alt_reads) / denominator if denominator else "NA",
        "alt_mapq_min": min(alt_mapq) if alt_mapq else "NA",
        "alt_mapq_max": max(alt_mapq) if alt_mapq else "NA",
        "alt_baseq_min": min(alt_baseq) if alt_baseq else "NA",
        "alt_baseq_max": max(alt_baseq) if alt_baseq else "NA",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--site", action="append", required=True, help="CHROM:POS:REF:ALT")
    parser.add_argument("--mapq", type=int, default=30)
    parser.add_argument("--baseq", type=int, default=30)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for specification in args.site:
        chrom, position, ref, alt = specification.split(":")
        rows.append(
            summarize_site(
                args.samtools,
                args.bam,
                args.sample,
                chrom,
                int(position),
                ref.upper(),
                alt.upper(),
                args.mapq,
                args.baseq,
            )
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
