#!/usr/bin/env python3
"""Compare two haploid mtDNA VCFs without treating missing calls as reference."""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


Call = tuple[str, int, str, str]


@dataclass(frozen=True)
class SpacerRecord:
    sample: str
    position: int
    ref: str
    alt: str
    gt: str
    dp: int
    idv: int
    imf: str
    alt_forward: int
    alt_reverse: int


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_info(value: str) -> dict[str, str]:
    output: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, entry = item.split("=", 1)
            output[key] = entry
    return output


def parse_vcf(path: Path, sample: str) -> tuple[set[Call], SpacerRecord]:
    calls: set[Call] = set()
    spacer: SpacerRecord | None = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, position_text, _identifier, ref, alt_text, _quality, _filter, info_text, fmt, value = fields[:10]
            position = int(position_text)
            format_values = dict(zip(fmt.split(":"), value.split(":")))
            gt = format_values.get("GT", ".").replace("|", "/")
            called_indices = {
                int(index)
                for index in gt.split("/")
                if index not in (".", "0")
            }
            alts = alt_text.split(",")
            for index in called_indices:
                if index < 1 or index > len(alts):
                    continue
                alt = alts[index - 1]
                if len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                    calls.add((chrom, position, ref, alt))

            spacer_index = alts.index("C") + 1 if "C" in alts else None
            if (
                position == 3106
                and ref == "CN"
                and spacer_index is not None
                and spacer_index in called_indices
            ):
                info = parse_info(info_text)
                dp4 = [int(item) for item in info["DP4"].split(",")]
                spacer = SpacerRecord(
                    sample=sample,
                    position=position,
                    ref=ref,
                    alt="C",
                    gt=gt,
                    dp=int(info["DP"]),
                    idv=int(info["IDV"]),
                    imf=info["IMF"],
                    alt_forward=dp4[2],
                    alt_reverse=dp4[3],
                )
    if spacer is None:
        raise ValueError(f"{sample}: expected called MT:3106 CN>C spacer record")
    return calls, spacer


def read_support(paths: list[Path]) -> dict[tuple[str, int], dict[str, str]]:
    output: dict[tuple[str, int], dict[str, str]] = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["sample"], int(row["position"]))
                if key in output:
                    raise ValueError(f"duplicate target-support row: {key}")
                output[key] = row
    for key in (("ORT15", 4215), ("ORT15", 7028), ("ORT16", 4215), ("ORT16", 7028)):
        if key not in output:
            raise ValueError(f"missing target-support row: {key}")
    return output


def validate_support(support: dict[tuple[str, int], dict[str, str]]) -> None:
    specifications = {
        ("ORT15", 4215): ("A", "G"),
        ("ORT15", 7028): ("C", "T"),
        ("ORT16", 4215): ("A", "G"),
        ("ORT16", 7028): ("C", "T"),
    }
    for key, (ref, alt) in specifications.items():
        row = support[key]
        if (row["chrom"], row["ref"], row["alt"]) != ("MT", ref, alt):
            raise ValueError(f"unexpected target-site alleles: {key}")
        if (row["mapq_min"], row["baseq_min"]) != ("30", "30"):
            raise ValueError(f"unexpected target-site quality filters: {key}")

    count_fields = (
        "qualifying_ref_reads",
        "qualifying_alt_reads",
        "qualifying_other_reads",
    )
    for key in (("ORT15", 4215), ("ORT15", 7028), ("ORT16", 4215)):
        if any(int(support[key][field]) != 0 for field in count_fields):
            raise ValueError(f"expected no qualifying coverage: {key}")

    row = support[("ORT16", 7028)]
    fields = (
        "qualifying_ref_reads",
        "qualifying_alt_reads",
        "distinct_alt_read_names",
        "alt_forward",
        "alt_reverse",
        "alt_fraction_ref_alt",
        "alt_mapq_min",
        "alt_mapq_max",
        "alt_baseq_min",
        "alt_baseq_max",
    )
    expected = ("0", "1", "1", "1", "0", "1.0", "37", "37", "40", "40")
    if tuple(row[field] for field in fields) != expected:
        raise ValueError("unexpected aggregate support for ORT16 C7028T")


def site_status(row: dict[str, str], call: Call, calls: set[Call]) -> str:
    covered = sum(
        int(row[field])
        for field in ("qualifying_ref_reads", "qualifying_alt_reads", "qualifying_other_reads")
    )
    if covered == 0:
        return "no_qualifying_coverage"
    if call in calls:
        return f"{call[2]}>{call[3]}_{int(row['qualifying_alt_reads'])}_read_call"
    return "covered_no_ALT_call"


def write_exact_calls(path: Path, calls: set[Call]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["CHROM", "POS", "REF", "ALT"])
        writer.writerows(sorted(calls, key=lambda row: (row[0], row[1], row[2], row[3])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ort15-vcf", type=Path, required=True)
    parser.add_argument("--ort16-vcf", type=Path, required=True)
    parser.add_argument("--target-support", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    ort15, spacer15 = parse_vcf(args.ort15_vcf, "ORT15")
    ort16, spacer16 = parse_vcf(args.ort16_vcf, "ORT16")
    support = read_support(args.target_support)
    validate_support(support)
    shared = ort15 & ort16
    only15 = ort15 - ort16
    only16 = ort16 - ort15
    if len(ort15) != 36 or len(ort16) != 37 or len(shared) != 36 or only15:
        raise ValueError(
            "corrected call-set invariant failed: expected ORT15/ORT16/shared=36/37/36"
        )
    if only16 != {("MT", 7028, "C", "T")}:
        raise ValueError(f"unexpected ORT16-only calls: {sorted(only16)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_exact_calls(args.output_dir / "ORT15.called_alt_snps.tsv", ort15)
    write_exact_calls(args.output_dir / "ORT16.called_alt_snps.tsv", ort16)
    write_exact_calls(args.output_dir / "ORT15_ORT16.shared_exact_snps.tsv", shared)

    pair_fields = [
        "record_type",
        "position",
        "ORT15",
        "ORT16",
        "ORT15_depth",
        "ORT16_depth",
        "interpretation",
    ]
    pair_rows = [
        {
            "record_type": "summary",
            "position": "NA",
            "ORT15": "36_called_ALT_SNPs",
            "ORT16": "37_called_ALT_SNPs",
            "ORT15_depth": "NA",
            "ORT16_depth": "NA",
            "interpretation": "36 exact shared POS_REF_ALT calls; ORT15 call set is a subset of ORT16",
        },
        {
            "record_type": "artifact",
            "position": 3106,
            "ORT15": "CN>C",
            "ORT16": "CN>C",
            "ORT15_depth": spacer15.dp,
            "ORT16_depth": spacer16.dp,
            "interpretation": "left-anchored deletion of artificial rCRS N spacer at 3107; exclude from biological SNP and discordance counts",
        },
        {
            "record_type": "site",
            "position": 4215,
            "ORT15": site_status(support[("ORT15", 4215)], ("MT", 4215, "A", "G"), ort15),
            "ORT16": site_status(support[("ORT16", 4215)], ("MT", 4215, "A", "G"), ort16),
            "ORT15_depth": sum(int(support[("ORT15", 4215)][field]) for field in ("qualifying_ref_reads", "qualifying_alt_reads", "qualifying_other_reads")),
            "ORT16_depth": sum(int(support[("ORT16", 4215)][field]) for field in ("qualifying_ref_reads", "qualifying_alt_reads", "qualifying_other_reads")),
            "interpretation": "old ORT15 A>G call does not survive; one truncated-read molecule was counted twice after nested-state merge",
        },
        {
            "record_type": "site",
            "position": 7028,
            "ORT15": site_status(support[("ORT15", 7028)], ("MT", 7028, "C", "T"), ort15),
            "ORT16": site_status(support[("ORT16", 7028)], ("MT", 7028, "C", "T"), ort16),
            "ORT15_depth": sum(int(support[("ORT15", 7028)][field]) for field in ("qualifying_ref_reads", "qualifying_alt_reads", "qualifying_other_reads")),
            "ORT16_depth": sum(int(support[("ORT16", 7028)][field]) for field in ("qualifying_ref_reads", "qualifying_alt_reads", "qualifying_other_reads")),
            "interpretation": "ORT16 call has one forward read (MAPQ37; BQ40; AF=1/1) and is not a confirmed pairwise difference",
        },
    ]
    with (args.output_dir / "mtdna_pair_concordance.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=pair_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(pair_rows)

    spacer_fields = [
        "sample",
        "position",
        "REF",
        "ALT",
        "GT",
        "DP",
        "IDV",
        "IMF",
        "ALT_forward",
        "ALT_reverse",
        "interpretation",
    ]
    with (args.output_dir / "mtdna_3106_spacer_audit.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=spacer_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in (spacer15, spacer16):
            writer.writerow(
                {
                    "sample": record.sample,
                    "position": record.position,
                    "REF": record.ref,
                    "ALT": record.alt,
                    "GT": record.gt,
                    "DP": record.dp,
                    "IDV": record.idv,
                    "IMF": record.imf,
                    "ALT_forward": record.alt_forward,
                    "ALT_reverse": record.alt_reverse,
                    "interpretation": "artificial rCRS N-spacer deletion; non-biological",
                }
            )


if __name__ == "__main__":
    main()
