#!/usr/bin/env python3
"""Plot fixed-comparator imputation concordance at manuscript display size.

Dosage r-squared uses column c5 of the GLIMPSE per-bin output; c4 is the
best-guess statistic. NRD is read from the independently harvested GCsV output,
checked against the manuscript summary, and plotted from a zero baseline.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RSQUARE = json.loads((HERE / "rsquare_all16.json").read_text())
GCSV = json.loads((HERE / "impute_gcsv.json").read_text())

# Colour-blind-aware colours plus distinct line styles and markers.
STRATEGIES = [
    ("Uncorrected", "raw", "#4D4D4D", ":", "x"),
    ("Trim-5", "trim5", "#E69F00", "-", "o"),
    ("Trim-10", "trim10", "#D55E00", "--", "s"),
    ("Rescale-5", "rescale5", "#0072B2", "-", "o"),
    ("Rescale-10", "rescale10", "#56B4E9", "--", "s"),
    ("Rescale-12", "rescaled", "#CC79A7", ":", "^"),
    ("bamRefine-5", "bamrefine5", "#009E73", "-", "D"),
    ("bamRefine-10", "bamrefine10", "#332288", "--", "P"),
]


def curve(sample: str, treatment: str) -> tuple[list[float], list[float]]:
    rows = [row for row in RSQUARE[f"{sample}_{treatment}"] if row["bin"] >= 1]
    return [row["meanMAF"] for row in rows], [row["c5"] for row in rows]


def weighted_r2(sample: str, treatment: str) -> tuple[float, float]:
    """Return count-weighted dosage and best-guess r-squared for bins 1--8."""
    rows = [
        row for row in RSQUARE[f"{sample}_{treatment}"]
        if 1 <= row["bin"] <= 8
    ]
    total = sum(int(row["n"]) for row in rows)
    if total == 0:
        raise ValueError(f"No MAF-bin 1--8 observations for {sample}/{treatment}")
    dosage = sum(int(row["n"]) * float(row["c5"]) for row in rows) / total
    best_guess = sum(int(row["n"]) * float(row["c4"]) for row in rows) / total
    return dosage, best_guess


def verify_r2_summary(
    summary: dict[tuple[str, str], dict[str, str]],
) -> int:
    """Check both published r-squared columns for all 16 comparisons."""
    checked = 0
    for sample in ("ORT15", "ORT16"):
        for _, treatment, _, _, _ in STRATEGIES:
            dosage, best_guess = weighted_r2(sample, treatment)
            table = summary[(sample, treatment)]
            expected = (
                float(table["dosage_r2_maf01"]),
                float(table["bestguess_r2_maf01"]),
            )
            for metric, value, target in zip(
                ("dosage", "best-guess"), (dosage, best_guess), expected,
            ):
                if round(value, 4) != round(target, 4):
                    raise AssertionError(
                        f"{metric} r2 mismatch {sample}/{treatment}: "
                        f"per-bin {value:.6f} vs table {target:.6f}"
                    )
                checked += 1
    return checked


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.07, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top")


def plot_curve(ax: plt.Axes, sample: str, title: str, label: str,
               show_ylabel: bool) -> None:
    for name, treatment, colour, line_style, marker in STRATEGIES:
        maf, dosage_r2 = curve(sample, treatment)
        ax.plot(maf, dosage_r2, color=colour, linestyle=line_style,
                marker=marker, markersize=4.2, linewidth=1.45, label=name)
    ax.set_xscale("log")
    ax.set_xlim(0.002, 0.6)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Minor allele frequency (MAF)")
    if show_ylabel:
        ax.set_ylabel("Dosage $r^2$")
    ax.set_title(title, fontsize=9.2, pad=7)
    ax.grid(axis="both", linestyle=":", linewidth=0.6, alpha=0.42)
    ax.tick_params(labelsize=8)
    panel_label(ax, label)


def plot_nrd(ax: plt.Axes, sample: str, title: str, label: str,
             summary: dict[tuple[str, str], dict[str, str]],
             show_ylabel: bool) -> int:
    names = [item[0] for item in STRATEGIES]
    keys = [item[1] for item in STRATEGIES]
    colours = [item[2] for item in STRATEGIES]
    values = [float(GCSV[sample][key][0]) for key in keys]

    # Compare the independent raw harvest with every manuscript table value.
    for key, value in zip(keys, values):
        table_value = float(summary[(sample, key)]["NRD_pct"])
        if round(value, 2) != table_value:
            raise AssertionError(
                f"NRD mismatch {sample}/{key}: GCsV {value} vs table {table_value}"
            )

    x = np.arange(len(names))
    bars = ax.bar(x, values, color=colours, edgecolor="#222222",
                  linewidth=0.45, width=0.72)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=42, ha="right", rotation_mode="anchor",
                       fontsize=7.2)
    # headroom for the rotated value labels: the tallest bar's label must clear
    # the top spine, so extend the axis above the data rather than clipping text
    ax.set_ylim(0, 27.0)
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    if show_ylabel:
        ax.set_ylabel("NRD (%)")
    ax.set_title(title, fontsize=9.2, pad=7)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.42)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.35,
                f"{value:.2f}", ha="center", va="bottom", fontsize=6.8,
                rotation=90)
    panel_label(ax, label)
    return len(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    records = list(csv.DictReader(args.summary.open(newline="", encoding="utf-8")))
    summary = {(row["sample"], row["treatment"]): row for row in records}
    if len(summary) != 16:
        raise ValueError(f"Expected 16 summary rows, found {len(summary)}")
    checked_r2 = verify_r2_summary(summary)

    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 6.4),
        gridspec_kw={"height_ratios": [1.15, 0.92]},
    )
    plot_curve(
        axes[0, 0], "ORT15",
        "ORT15: non-UDG 0.27$\\times$ vs full-UDG 0.90$\\times$", "a", True,
    )
    plot_curve(
        axes[0, 1], "ORT16",
        "ORT16: non-UDG 0.35$\\times$ vs full-UDG 1.46$\\times$", "b", False,
    )
    checked = plot_nrd(axes[1, 0], "ORT15", "ORT15", "c", summary, True)
    checked += plot_nrd(axes[1, 1], "ORT16", "ORT16", "d", summary, False)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               fontsize=8, handlelength=2.6, columnspacing=1.25,
               bbox_to_anchor=(0.5, 0.006))
    fig.subplots_adjust(left=0.10, right=0.985, top=0.965, bottom=0.185,
                        hspace=0.42, wspace=0.22)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        args.output,
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None,
                  "Creator": "ORT PMD deterministic Figure 3 generator"},
    )
    plt.close(fig)

    rows = []
    for sample in ("ORT15", "ORT16"):
        for _, treatment, _, _, _ in STRATEGIES:
            for row in RSQUARE[f"{sample}_{treatment}"]:
                if row["bin"] >= 1:
                    rows.append((sample, treatment, row["bin"],
                                 round(row["meanMAF"], 6), round(row["c5"], 6),
                                 round(row["c4"], 6)))
    plot_data = args.output.parent / "figure_03_plotdata_perbin.tsv"
    with plot_data.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "sample", "treatment", "maf_bin", "mean_MAF",
            "dosage_r2", "bestguess_r2",
        ])
        writer.writerows(rows)

    print("=== Figure 3 fixed-comparator verification ===")
    print(f"weighted dosage/best-guess r2 == table CSV (4 dp) for all {checked_r2}: PASS")
    print(f"NRD GCsV harvest == table CSV (2 dp) for all {checked}: PASS")
    print(f"per-bin plot-data rows: {len(rows)} (expected 128): PASS")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
