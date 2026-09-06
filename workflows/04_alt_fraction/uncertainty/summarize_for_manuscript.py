#!/usr/bin/env python3
"""Build manuscript-facing DeltaALT summaries from validated jackknife outputs.

This script does not read the primary per-site TSVs and does not alter any
manuscript.  It reduces the validated state/contrast tables to the rows needed
for Figure 2, Supplementary Table S2, and the Results text, and derives paired
rescaling-window contrasts from identically omitted genomic blocks.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: (f"{value:.12g}" if isinstance(value, float) else value)
                for key, value in row.items()
            })


def jackknife(full: float, deleted: list[float], critical: float) -> dict[str, float]:
    n = len(deleted)
    mean_deleted = math.fsum(deleted) / n
    se = math.sqrt((n - 1) / n * math.fsum((x - mean_deleted) ** 2 for x in deleted))
    bias_corrected = n * full - (n - 1) * mean_deleted
    return {
        "full_estimate_pp": full,
        "jackknife_bias_corrected_pp": bias_corrected,
        "jackknife_se_pp": se,
        "ci95_lo_pp": bias_corrected - critical * se,
        "ci95_hi_pp": bias_corrected + critical * se,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("primary_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    states = read_rows(args.primary_dir / "jackknife_state_summary.csv")
    state_delete = read_rows(args.primary_dir / "jackknife_state_delete_one.csv")
    contrasts = read_rows(args.primary_dir / "jackknife_contrast_summary.csv")

    loco_states = [row for row in states if row["scheme"] == "loco"]
    if len(loco_states) != 26:
        raise AssertionError(f"Expected 26 LOCO state rows, observed {len(loco_states)}")
    state_out: list[dict[str, object]] = []
    for row in loco_states:
        point = float(row["full_estimate_pp"])
        lo = float(row["ci95_lo_pp"])
        hi = float(row["ci95_hi_pp"])
        state_out.append({
            "individual": row["individual"],
            "library": row["library"],
            "treatment": row["treatment"],
            "deltaalt_point_pp": point,
            "ci95_lo_pp": lo,
            "ci95_hi_pp": hi,
            "display_2dp": f"{point:+.2f} ({lo:+.2f} to {hi:+.2f})",
            "interval_excludes_zero": row["ci_excludes_zero"],
            "uncertainty_unit": "leave-one-autosome-out_95pct_CI",
        })
    write_rows(args.output_dir / "manuscript_state_ci_loco.csv", state_out)

    selected_ids = {
        "raw_nu_minus_raw_fu",
        "nu_bamrefine5_minus_rescale5",
        "rescale5_shift_nu_minus_fu",
        "rescale12_shift_nu_minus_fu",
    }
    selected = [
        row for row in contrasts
        if row["scheme"] == "loco"
        and (row["contrast_id"] in selected_ids or "_minus_raw" in row["contrast_id"])
    ]
    contrast_out: list[dict[str, object]] = []
    for row in selected:
        point = float(row["full_estimate_pp"])
        lo = float(row["ci95_lo_pp"])
        hi = float(row["ci95_hi_pp"])
        contrast_out.append({
            "individual": row["individual"],
            "contrast_id": row["contrast_id"],
            "description": row["description"],
            "estimate_pp": point,
            "ci95_lo_pp": lo,
            "ci95_hi_pp": hi,
            "display_2dp": f"{point:+.2f} ({lo:+.2f} to {hi:+.2f})",
            "interval_excludes_zero": row["ci_excludes_zero"],
        })
    write_rows(args.output_dir / "manuscript_primary_contrasts_loco.csv", contrast_out)

    summary_index = {
        (row["scheme"], row["individual"], row["library"], row["treatment"]): row
        for row in states
    }
    delete_index: dict[tuple[str, str, str, str], dict[tuple[str, str], float]] = defaultdict(dict)
    for row in state_delete:
        key = (row["scheme"], row["individual"], row["library"], row["treatment"])
        block = (row["omitted_chrom"], row["omitted_block_index_0based"])
        delete_index[key][block] = float(row["delete_one_deltaalt_pp"])

    window_out: list[dict[str, object]] = []
    for scheme in ("loco", "5mb", "10mb"):
        for individual in ("ORT15", "ORT16"):
            for library in ("fu", "nu"):
                longer_windows = ["rescale12"]
                if library == "nu":
                    longer_windows.append("rescale10")
                for long_name in longer_windows:
                    short_key = (scheme, individual, library, "rescale5")
                    long_key = (scheme, individual, library, long_name)
                    a = summary_index[long_key]
                    b = summary_index[short_key]
                    blocks = sorted(set(delete_index[long_key]) | set(delete_index[short_key]))
                    if set(delete_index[long_key]) != set(delete_index[short_key]):
                        raise AssertionError("Paired rescaling states have different block universes")
                    deleted = [delete_index[long_key][block] - delete_index[short_key][block] for block in blocks]
                    full = float(a["full_estimate_pp"]) - float(b["full_estimate_pp"])
                    values = jackknife(full, deleted, float(a["t_critical_975"]))
                    window_out.append({
                        "scheme": scheme,
                        "individual": individual,
                        "library": library,
                        "contrast_id": f"{library}_{long_name}_minus_rescale5",
                        "direction": f"{long_name} minus rescale5",
                        "n_blocks": len(blocks),
                        **values,
                        "ci_excludes_zero": not (
                            values["ci95_lo_pp"] <= 0 <= values["ci95_hi_pp"]
                        ),
                    })
    write_rows(args.output_dir / "paired_rescaling_window_contrasts.csv", window_out)


if __name__ == "__main__":
    main()
