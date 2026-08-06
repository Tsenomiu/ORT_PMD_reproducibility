#!/usr/bin/env python3
"""Summarize native ancIBD IBD1 segments using strict length/density cutoffs."""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

PATTERN = re.compile(r"^(ORT15|ORT16)_(fu|nu)_(raw|trim5|rescale5|bamrefine5)$")
THRESHOLDS = (8, 12, 16, 20)
DENOMINATOR_CM = 3538.8511


def meta(iid1: str, iid2: str) -> dict[str, str]:
    a, b = PATTERN.fullmatch(iid1), PATTERN.fullmatch(iid2)
    if not a or not b or a.group(1) != "ORT15" or b.group(1) != "ORT16":
        raise ValueError(f"Unexpected oriented pair: {iid1}, {iid2}")
    same = a.group(3) == b.group(3); mixed = a.group(2) != b.group(2)
    non_udg = a.group(3) if a.group(2) == "nu" and b.group(2) == "fu" else b.group(3) if a.group(2) == "fu" and b.group(2) == "nu" else ""
    return {"iid1": iid1, "iid2": iid2, "individual1": "ORT15", "individual2": "ORT16",
            "library1": a.group(2), "library2": b.group(2), "treatment1": a.group(3), "treatment2": b.group(3),
            "library_direction": f"ORT15_{a.group(2)}_x_ORT16_{b.group(2)}", "same_treatment": "yes" if same else "no",
            "mixed_library": "yes" if mixed else "no", "same_treatment_method": a.group(3) if same else "", "non_udg_treatment": non_udg}


def write(path: Path, rows: list[dict]) -> None:
    if not rows: raise ValueError(f"No rows for {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("segments", type=Path)
    parser.add_argument("pair_design", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.pair_design.open(newline="") as handle:
        pairs = [(r["iid1"], r["iid2"]) for r in csv.DictReader(handle, delimiter="\t")]
    if len(pairs) != 64 or len(set(pairs)) != 64: raise ValueError("Expected 64 unique pairs")
    segments = []; by_pair = defaultdict(list)
    with args.segments.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["iid1"], row["iid2"])
            if key not in set(pairs) and key[::-1] in set(pairs): key = key[::-1]
            length_cm = 100 * float(row["lengthM"]); marker_steps = int(float(row["length"])); density = marker_steps / length_cm
            item = {**meta(*key), "chromosome": int(row["ch"]), "start_cM": f"{100*float(row['StartM']):.8f}",
                    "end_cM": f"{100*float(row['EndM']):.8f}", "length_cM": f"{length_cm:.8f}",
                    "start_bp": int(row["StartBP"]), "end_bp": int(row["EndBP"]), "marker_steps": marker_steps,
                    "snp_per_cM": f"{density:.8f}", "segment_type": "IBD1"}
            segments.append(item); by_pair[key].append((length_cm, density))
    summaries = []
    for key in pairs:
        item = {**meta(*key), "callable_map_denominator_cM": f"{DENOMINATOR_CM:.6f}"}
        for label, density_cut in (("native", 0), ("density220", 220)):
            for threshold in THRESHOLDS:
                selected = [length for length, density in by_pair[key] if length > threshold and density > density_cut]
                total = sum(selected)
                item[f"{label}_sum_gt{threshold}_cM"] = f"{total:.6f}"
                item[f"{label}_n_gt{threshold}"] = len(selected)
                item[f"{label}_recovery_gt{threshold}_pct"] = f"{100*total/DENOMINATOR_CM:.6f}"
        summaries.append(item)
    write(args.output_dir / "CROSS64_PAIR_SUMMARY.tsv", summaries)
    write(args.output_dir / "CROSS64_SEGMENTS_NATIVE.tsv", segments)
    write(args.output_dir / "CROSS64_SEGMENTS_DENSITY220.tsv", [r for r in segments if float(r["length_cM"]) > 8 and float(r["snp_per_cM"]) > 220])
    mixed = [r for r in summaries if r["same_treatment"] == "yes" and r["mixed_library"] == "yes"]
    asymmetric = [r for r in summaries if (r["library1"] == "fu" and r["treatment1"] == "raw" and r["library2"] == "nu") or (r["library1"] == "nu" and r["library2"] == "fu" and r["treatment2"] == "raw")]
    if len(mixed) != 8 or len(asymmetric) != 8: raise ValueError("Paper-design subset must contain 8 rows each")
    write(args.output_dir / "CROSS64_SAME_TREATMENT_MIXED.tsv", mixed)
    write(args.output_dir / "CROSS64_ASYMMETRIC_FULLUDG_RAW.tsv", asymmetric)


if __name__ == "__main__":
    main()

