#!/usr/bin/env python3
"""Summarize ALT fractions (DP>=3) and retention (DP>0) without mixing denominators."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def site_class(ref: str, alt: str) -> str | None:
    ref, alt = ref.upper(), alt.upper()
    if (ref, alt) in {("C", "T"), ("G", "A")}:
        return "damage_transition"
    if ref in "ACGT" and alt in "ACGT" and ref != alt:
        if {ref, alt} not in ({"A", "G"}, {"C", "T"}):
            return "transversion"
    return None


def summarize(path: Path) -> dict[str, float | int]:
    sums = defaultdict(float)
    n = defaultdict(int)
    covered = 0
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dp = int(row["DP"] or 0)
            ad_ref = int((row["AD_REF"] or "0").split(",")[0])
            ad_alt = int((row["AD_ALT"] or "0").split(",")[0])
            depth = ad_ref + ad_alt
            if depth > 0:
                covered += 1
            cls = site_class(row["REF"], row["ALT"])
            if cls and depth >= 3:
                sums[cls] += ad_alt / depth
                n[cls] += 1
    return {
        "covered_sites_any_depth": covered,
        "damage_transition_n_dp3": n["damage_transition"],
        "damage_transition_mean_pct_dp3": 100 * sums["damage_transition"] / n["damage_transition"] if n["damage_transition"] else float("nan"),
        "transversion_n_dp3": n["transversion"],
        "transversion_mean_pct_dp3": 100 * sums["transversion"] / n["transversion"] if n["transversion"] else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    values = {}
    for item in manifest:
        key = (item["sample"], item["library"], item["treatment"])
        values[key] = summarize(Path(item["counts_tsv"]))
    for item in manifest:
        key = (item["sample"], item["library"], item["treatment"])
        value = values[key]
        raw = values[(item["sample"], item["library"], "raw")]["covered_sites_any_depth"]
        rows.append({
            **{k: item[k] for k in ("sample", "library", "treatment")},
            **value,
            "retention_pct_all_covered": 100 * int(value["covered_sites_any_depth"]) / int(raw),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()

