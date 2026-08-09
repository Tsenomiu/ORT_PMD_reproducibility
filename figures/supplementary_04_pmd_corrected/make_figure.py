#!/usr/bin/env python3
"""Generate corrected PMD profiles from included mapDamage tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"


def profile(path: Path, end: str, change: str, reference: str):
    points = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["End"] != end or row["Std"] != "+":
                continue
            pos = int(row["Pos"])
            if not 1 <= pos <= 25:
                continue
            denominator = float(row[reference])
            value = float(row[change]) / denominator if denominator else float("nan")
            points.append((pos if end == "5p" else -pos, value))
    return tuple(zip(*sorted(points)))


def input_path(individual: str, library: str, treatment: str) -> Path:
    prefix = f"{individual}_{library}"
    if treatment == "rescale5":
        # Rescaling changes base qualities, not mismatch counts; the mapDamage
        # profile is therefore the same as the uncorrected aggregate table.
        treatment = "uncorrected"
    return INPUTS / f"{prefix}_{treatment}.txt"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    treatments = (
        ("Trim-5", "trim5", "#1F77B4", "-"),
        ("Rescale-5", "rescale5", "#2CA02C", "-"),
        ("bamRefine-5", "bamrefine5", "#D62728", "-"),
        ("Uncorrected", "uncorrected", "#4D4D4D", "--"),
    )
    rows = (
        ("ORT15", "fullUDG", "ORT15 full-UDG (merged)"),
        ("ORT15", "nonUDG", "ORT15 non-UDG (post-dedup)"),
        ("ORT16", "fullUDG", "ORT16 full-UDG (merged)"),
        ("ORT16", "nonUDG", "ORT16 non-UDG (post-dedup)"),
    )

    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )
    fig, axes = plt.subplots(4, 2, figsize=(13, 14.2), sharey=True)
    for row_index, (individual, library, row_title) in enumerate(rows):
        for column_index, (end, change, reference, title, xlabel, ylabel) in enumerate(
            (
                ("5p", "C>T", "C", "5' C->T", "Position from 5' end (bp)", "C->T frequency"),
                ("3p", "G>A", "G", "3' G->A", "Position from 3' end (bp)", "G->A frequency"),
            )
        ):
            axis = axes[row_index, column_index]
            for label, key, color, linestyle in treatments:
                x, y = profile(input_path(individual, library, key), end, change, reference)
                axis.plot(x, y, label=label, color=color, linestyle=linestyle, linewidth=1.4)
            if column_index == 0:
                axis.set_title(f"{row_title}             {title}", loc="left", fontweight="bold")
            else:
                axis.set_title(title, fontweight="bold")
            axis.set_ylabel(ylabel)
            if row_index == 3:
                axis.set_xlabel(xlabel)
            axis.set_ylim(0, 0.30)
            axis.set_yticks((0, 0.1, 0.2, 0.3))
            axis.set_xticks((1, 5, 9, 13, 17, 21, 25) if end == "5p" else (-25, -21, -17, -13, -9, -5, -1))
            axis.grid(True, linestyle=":", color="#B8B8B8", alpha=0.7)
            if row_index == 0 and column_index == 0:
                axis.legend(loc="upper right", frameon=True, fontsize=8)

    fig.suptitle("mapDamage misincorporation profiles after PMD correction (collapsed BAMs)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.965), h_pad=1.4)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
