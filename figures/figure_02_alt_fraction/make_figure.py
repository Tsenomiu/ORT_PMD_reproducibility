#!/usr/bin/env python3
"""Plot DeltaALT, LOCO confidence intervals, and covered-site retention."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from PIL import Image


COLORS = {"ORT15": "#2B6CB0", "ORT16": "#C43D3D"}
MARKERS = {"ORT15": "o", "ORT16": "^"}
LABELS = {
    "raw": "No correction",
    "trim5": "Trim-5",
    "trim10": "Trim-10",
    "rescale5": "Rescale-5",
    "rescale10": "Rescale-10",
    "rescale12": "Rescale-12",
    "bamrefine5": "bamRefine-5",
    "bamrefine10": "bamRefine-10",
}
ORDER = {
    "fu": ["raw", "trim5", "rescale5", "rescale12", "bamrefine5"],
    "nu": [
        "raw", "trim5", "trim10", "rescale5", "rescale10",
        "rescale12", "bamrefine5", "bamrefine10",
    ],
}


def expected_keys() -> set[tuple[str, str, str]]:
    return {
        (individual, library, treatment)
        for individual in COLORS
        for library, treatments in ORDER.items()
        for treatment in treatments
    }


def read_rows(path: Path) -> dict[tuple[str, str, str], dict[str, float]]:
    rows: dict[tuple[str, str, str], dict[str, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["individual"], row["library"], row["treatment"])
            if key not in expected_keys():
                raise ValueError(f"Unexpected state in point-estimate table: {key}")
            if key in rows:
                raise ValueError(f"Duplicate state in point-estimate table: {key}")
            rows[key] = {
                "delta": float(row["ts_minus_tv_pp"]),
                "retention": float(row["retention_pct"]),
            }
    if set(rows) != expected_keys():
        raise ValueError("Point-estimate table does not contain the exact 26 expected states")
    return rows


def read_intervals(path: Path) -> dict[tuple[str, str, str], tuple[float, float, float]]:
    rows: dict[tuple[str, str, str], tuple[float, float, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["individual"], row["library"], row["treatment"])
            if key not in expected_keys():
                raise ValueError(f"Unexpected state in interval table: {key}")
            if key in rows:
                raise ValueError(f"Duplicate state in interval table: {key}")
            if row.get("uncertainty_unit") != "leave-one-autosome-out_95pct_CI":
                raise ValueError(f"Interval is not a LOCO 95% CI for {key}")
            rows[key] = (
                float(row["deltaalt_point_pp"]),
                float(row["ci95_lo_pp"]),
                float(row["ci95_hi_pp"]),
            )
    if set(rows) != expected_keys():
        raise ValueError("Interval table does not contain the exact 26 expected states")
    return rows


def draw_panel(ax: plt.Axes, rows, intervals, library: str, title: str, panel: str) -> None:
    order = ORDER[library]
    x = np.arange(len(order), dtype=float)
    bar_width = 0.32
    ax_ret = ax.twinx()

    for individual, offset in (("ORT15", -bar_width / 2), ("ORT16", bar_width / 2)):
        color = COLORS[individual]
        delta = [rows[(individual, library, treatment)]["delta"] for treatment in order]
        ci = [intervals[(individual, library, treatment)] for treatment in order]
        for treatment, point, interval in zip(order, delta, ci):
            # Figure point estimates are the canonical six-decimal values;
            # interval input retains full precision from the same estimator.
            if abs(point - interval[0]) > 5.1e-7:
                raise ValueError(
                    f"Point estimate mismatch for {individual}/{library}/{treatment}: "
                    f"{point} versus {interval[0]}"
                )
        lower = [point - interval[1] for point, interval in zip(delta, ci)]
        upper = [interval[2] - point for point, interval in zip(delta, ci)]
        if min(lower + upper) < 0:
            raise ValueError("A point estimate lies outside its reported confidence interval")
        retention = [rows[(individual, library, treatment)]["retention"] for treatment in order]
        ax_ret.bar(
            x + offset,
            retention,
            width=bar_width,
            color=color,
            alpha=0.13,
            edgecolor=color,
            linewidth=0.6,
            zorder=1,
        )
        ax.errorbar(
            x,
            delta,
            yerr=np.asarray([lower, upper]),
            color=color,
            fmt=MARKERS[individual] + "-",
            markersize=6,
            linewidth=1.8,
            elinewidth=0.85,
            capsize=2.3,
            capthick=0.85,
            zorder=4,
        )

    ax.axhline(0, color="#444444", linewidth=1.0, linestyle="--", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[t] for t in order], rotation=30, ha="right")
    ax.set_ylim(-2.0, 1.1)
    ax.set_yticks(np.arange(-2.0, 1.1, 0.5))
    ax.grid(axis="y", color="#D9D9D9", linestyle=":", linewidth=0.7, zorder=0)
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=7)
    ax.text(-0.12, 1.08, panel, transform=ax.transAxes, fontsize=12, fontweight="bold")
    ax.tick_params(axis="both", labelsize=8)
    ax.spines["top"].set_visible(False)

    ax_ret.set_ylim(0, 112)
    ax_ret.set_yticks([0, 20, 40, 60, 80, 100])
    ax_ret.tick_params(axis="y", labelsize=7, colors="#666666")
    ax_ret.spines["top"].set_visible(False)
    ax_ret.spines["right"].set_color("#999999")
    if library == "nu":
        ax_ret.set_ylabel("Covered SNP sites retained (%)", color="#666666", fontsize=8.5)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("intervals", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = read_rows(args.data)
    intervals = read_intervals(args.intervals)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.55), sharey=True)
    draw_panel(axes[0], rows, intervals, "fu", "Full-UDG libraries", "a")
    draw_panel(axes[1], rows, intervals, "nu", "Non-UDG libraries", "b")
    axes[0].set_ylabel(
        "Transition $-$ transversion ALT fraction\n(percentage points)",
        fontsize=8.7,
    )

    legend = [
        Line2D([0], [0], color=COLORS["ORT15"], marker=MARKERS["ORT15"],
               linewidth=1.8, markersize=6, label="ORT15"),
        Line2D([0], [0], color=COLORS["ORT16"], marker=MARKERS["ORT16"],
               linewidth=1.8, markersize=6, label="ORT16"),
        Patch(facecolor="#777777", edgecolor="#777777", alpha=0.13,
              label="Pale bars: covered-site retention"),
        Line2D([0], [0], color="#444444", linestyle="--", linewidth=1,
               label="Equal transition and transversion fractions"),
    ]
    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=False,
        fontsize=8,
        columnspacing=1.4,
    )
    fig.subplots_adjust(left=0.105, right=0.93, top=0.91, bottom=0.27, wspace=0.18)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    word_png = args.output.with_name(args.output.stem + "_word.png")
    fig.savefig(word_png, dpi=200, bbox_inches="tight", facecolor="white")
    with Image.open(word_png) as raster:
        rgba = raster.convert("RGBA")
        rgb = Image.new("RGB", rgba.size, "white")
        rgb.paste(rgba, mask=rgba.getchannel("A"))
        rgb.save(word_png, dpi=(200, 200), optimize=True)
    print(args.output)
    print(word_png)


if __name__ == "__main__":
    main()
