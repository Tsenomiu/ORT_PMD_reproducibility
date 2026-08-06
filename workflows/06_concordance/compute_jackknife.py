#!/usr/bin/env python3
"""Paired 22-autosome delete-one jackknife for bamRefine-5 minus Rescale-5 NRD."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


REFERENCE = {
    ("ORT15", "bamrefine5"): 18.30,
    ("ORT15", "rescale5"): 18.32,
    ("ORT16", "bamrefine5"): 14.88,
    ("ORT16", "rescale5"): 14.85,
}


def counts(fields: list[str]) -> tuple[int, int]:
    # GCsV zero-based fields: 8 RA_ok, 9 AA_ok, 10--12 three error classes.
    ra_ok, aa_ok = int(fields[8]), int(fields[9])
    errors = sum(int(fields[i]) for i in (10, 11, 12))
    return errors, ra_ok + aa_ok + errors


def nrd(rows: dict[int, tuple[int, int]], chroms: list[int]) -> float:
    numerator = sum(rows[c][0] for c in chroms)
    denominator = sum(rows[c][1] for c in chroms)
    return 100.0 * numerator / denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gcsv_rows", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    data: dict[tuple[str, str], dict[int, tuple[int, int]]] = {}
    with args.gcsv_rows.open() as handle:
        for line in handle:
            if "|" not in line:
                continue
            label, record = line.split("|", 1)
            individual, treatment, chrom_text = label.split()
            fields = record.split()
            if fields and fields[0] == "GCsV":
                data.setdefault((individual, treatment), {})[int(chrom_text)] = counts(fields)
    expected = set(REFERENCE)
    if set(data) != expected or any(set(v) != set(range(1, 23)) for v in data.values()):
        raise ValueError("Expected exactly four datasets with chromosomes 1--22")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "perchr_nrd_counts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["individual", "treatment", "chr", "nrd_numerator", "nrd_denominator", "nrd_pct"])
        for key in sorted(data):
            for chrom in range(1, 23):
                num, den = data[key][chrom]
                writer.writerow([*key, chrom, num, den, f"{100*num/den:.4f}"])
    chroms = list(range(1, 23)); tcrit = 2.080
    output = []
    for key, reference in REFERENCE.items():
        observed = nrd(data[key], chroms)
        if abs(observed - reference) >= 0.01:
            raise ValueError(f"{key} count sum {observed:.4f} differs from table {reference:.2f}")
    for individual in ("ORT15", "ORT16"):
        bam = data[(individual, "bamrefine5")]
        res = data[(individual, "rescale5")]
        full = nrd(bam, chroms) - nrd(res, chroms)
        deleted = [nrd(bam, [x for x in chroms if x != c]) - nrd(res, [x for x in chroms if x != c]) for c in chroms]
        mean_deleted = sum(deleted) / 22
        se = math.sqrt(21 / 22 * sum((x - mean_deleted) ** 2 for x in deleted))
        estimate = 22 * full - 21 * mean_deleted
        lo, hi = estimate - tcrit * se, estimate + tcrit * se
        output.append([individual, f"{nrd(bam, chroms):.4f}", f"{nrd(res, chroms):.4f}", f"{full:+.4f}", f"{estimate:+.4f}", f"{se:.4f}", f"{lo:+.4f}", f"{hi:+.4f}", "no" if lo < 0 < hi else "yes"])
    with (args.output_dir / "jackknife_deltaNRD.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["individual", "NRD_bamRefine5_pct", "NRD_Rescale5_pct", "deltaNRD_pp", "jackknife_est_pp", "jackknife_SE_pp", "CI95_lo_pp", "CI95_hi_pp", "CI_excludes_zero"])
        writer.writerows(output)


if __name__ == "__main__":
    main()
