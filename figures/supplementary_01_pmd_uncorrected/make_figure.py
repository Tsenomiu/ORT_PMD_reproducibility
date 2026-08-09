#!/usr/bin/env python3
"""Generate the uncorrected PMD-profile figure from mapDamage tables."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = {
        ("ORT15", "full-UDG"): INPUTS / "ORT15_fullUDG.txt",
        ("ORT15", "non-UDG"): INPUTS / "ORT15_nonUDG.txt",
        ("ORT16", "full-UDG"): INPUTS / "ORT16_fullUDG.txt",
        ("ORT16", "non-UDG"): INPUTS / "ORT16_nonUDG.txt",
    }
    colors = {"full-UDG": "#08306B", "non-UDG": "#D62728"}

    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharey=True)
    for row_index, individual in enumerate(("ORT15", "ORT16")):
        for column_index, (end, change, reference, title, xlabel, ylabel) in enumerate(
            (
                ("5p", "C>T", "C", "5' C->T", "Position from 5' end (bp)", "C->T frequency"),
                ("3p", "G>A", "G", "3' G->A", "Position from 3' end (bp)", "G->A frequency"),
            )
        ):
            axis = axes[row_index, column_index]
            for library in ("full-UDG", "non-UDG"):
                x, y = profile(paths[(individual, library)], end, change, reference)
                axis.plot(x, y, color=colors[library], linewidth=2.0, label=f"{library} merged")
            axis.set_title(f"{individual}  {title}", fontweight="bold")
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)
            axis.set_ylim(0, 0.30)
            axis.set_yticks((0, 0.1, 0.2, 0.3))
            axis.set_xticks((1, 5, 9, 13, 17, 21, 25) if end == "5p" else (-25, -21, -17, -13, -9, -5, -1))
            axis.grid(True, linestyle=":", color="#B8B8B8", alpha=0.7)
            if column_index == 0:
                axis.legend(loc="upper right", frameon=True, fontsize=8)

    fig.suptitle("mapDamage misincorporation profiles", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
