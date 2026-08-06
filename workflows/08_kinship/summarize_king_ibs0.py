#!/usr/bin/env python3
"""KING-robust kinship and phasing-independent IBS0 from a two-sample query TSV."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def hom(gt: str) -> int | None:
    gt = gt.replace("|", "/")
    return 0 if gt == "0/0" else 1 if gt == "1/1" else None


def het(gt: str) -> bool:
    return gt.replace("|", "/") in {"0/1", "1/0"}


def max_gp(text: str) -> float:
    return max(float(x) for x in text.split(",") if x != ".")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-max-gp", type=float)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    n = ibs0 = het1 = het2 = both_het = 0; unrelated_sum = 0.0
    with args.input.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.split()
            if len(fields) != 7:
                raise ValueError(f"Expected CHROM POS RAF GT GP GT GP at line {line_number}")
            _, _, p_text, gt1, gp1, gt2, gp2 = fields
            if args.min_max_gp is not None and (max_gp(gp1) < args.min_max_gp or max_gp(gp2) < args.min_max_gp):
                continue
            p = float(p_text.split(",")[0]); q = 1 - p
            h1, h2 = hom(gt1), hom(gt2)
            ibs0 += int(h1 is not None and h2 is not None and h1 != h2)
            ahet, bhet = het(gt1), het(gt2)
            het1 += ahet; het2 += bhet; both_het += ahet and bhet
            unrelated_sum += 2 * p * p * q * q; n += 1
    if not n or not (het1 + het2): raise ValueError("No informative sites")
    unrelated = unrelated_sum / n; full_sib = unrelated / 4
    observed = ibs0 / n; king = (both_het - 2 * ibs0) / (het1 + het2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sites", "IBS0", "KING_robust", "observed_IBS0_per_site", "unrelated_expected", "full_sibling_expected", "observed_over_full_sibling", "min_max_GP"])
        writer.writerow([n, ibs0, f"{king:.8f}", f"{observed:.8f}", f"{unrelated:.8f}", f"{full_sib:.8f}", f"{observed/full_sib:.8f}", args.min_max_gp or "none"])


if __name__ == "__main__":
    main()

