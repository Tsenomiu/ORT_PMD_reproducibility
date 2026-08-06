#!/usr/bin/env python3
"""Compute published weighted r2 (MAF bins 1--8) and all-site GLIMPSE NRD."""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def lines(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        yield from handle


def weighted_r2(path: Path) -> tuple[float, float, int]:
    dosage_num = best_num = 0.0
    total = 0
    for raw in lines(path):
        fields = raw.split()
        if not fields or fields[0].startswith("#") or not fields[0].isdigit():
            continue
        bin_number, count = int(fields[0]), int(float(fields[1]))
        if 1 <= bin_number <= 8 and count:
            # GLIMPSE2 rsquare.grp columns: bin, n, mean_AF, best-guess r2,
            # dosage r2. The final table must retain this ordering (for example,
            # ORT15 raw dosage 0.4354 > best guess 0.4175).
            best_num += count * float(fields[3])
            dosage_num += count * float(fields[4])
            total += count
    if not total:
        raise ValueError(f"No MAF-bin 1--8 observations in {path}")
    return dosage_num / total, best_num / total, total


def nrd(path: Path) -> float:
    for raw in lines(path):
        fields = raw.split()
        if fields and fields[0] == "GCsV":
            ra_ok, aa_ok = int(fields[8]), int(fields[9])
            errors = sum(int(fields[i]) for i in (10, 11, 12))
            return 100.0 * errors / (ra_ok + aa_ok + errors)
    raise ValueError(f"No GCsV row in {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    with args.manifest.open(newline="") as handle:
        for item in csv.DictReader(handle, delimiter="\t"):
            dosage, best, count = weighted_r2(Path(item["rsquare_grp"]))
            rows.append({
                "individual": item["individual"], "treatment": item["treatment"],
                "dosage_r2_maf01": f"{dosage:.6f}", "bestguess_r2_maf01": f"{best:.6f}",
                "n_sites_maf01": count,
                "NRD_pct_all_sites": f"{nrd(Path(item['error_spl'])):.4f}",
            })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
