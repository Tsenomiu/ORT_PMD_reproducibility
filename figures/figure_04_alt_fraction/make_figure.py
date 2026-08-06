#!/usr/bin/env python3
"""Plot alternative-allele fractions and retained sites by PMD treatment.

Two-row x two-column figure.
  Top row    : full-UDG (ORT15, ORT16)
  Bottom row : non-UDG  (ORT15, ORT16)
  y-axis     : mean alternative-allele fraction (%)
  x-axis     : PMD-handling strategy (ordered by aggressiveness)
  Filled circles  -> C->T / G->A transitions (PMD-prone)
  Open squares    -> transversions (PMD-robust baseline; only known for Untrimmed)
  Horizontal dashed line = untrimmed transversion baseline.

Numerical values are read from a tab-separated summary table.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Strategy ordering (left -> right)
ORDER_FU = ["Untrimmed", "Trim-5", "Rescale-5", "Rescale-12", "bamRefine-5"]
ORDER_NU = ["Untrimmed", "Trim-5", "Trim-10",
            "Rescale-5", "Rescale-10", "Rescale-12",
            "bamRefine-5", "bamRefine-10"]

TREATMENT_LABELS = {
    "raw": "Untrimmed",
    "trim5": "Trim-5",
    "trim10": "Trim-10",
    "rescale5": "Rescale-5",
    "rescale10": "Rescale-10",
    "rescale12": "Rescale-12",
    "bamrefine5": "bamRefine-5",
    "bamrefine10": "bamRefine-10",
}

COLOR = {"ORT15": "#1f77b4", "ORT16": "#d62728"}


def load_data(path: Path) -> list[dict[str, object]]:
    data = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            treatment = row["treatment"]
            if treatment not in TREATMENT_LABELS:
                raise ValueError(f"Unknown treatment label: {treatment}")
            data.append({
                "ind": row["individual"],
                "lib": row["library"],
                "strategy": TREATMENT_LABELS[treatment],
                "ts": float(row["damage_transition_mean_pct_dp3"]),
                "tv": (float("nan") if row["transversion_mean_pct_dp3"] == "NA"
                       else float(row["transversion_mean_pct_dp3"])),
                "sites_M": float(row["covered_sites_any_depth_millions"]),
                "retain": float(row["retention_pct_all_covered"]),
            })
    return data


def panel(ax, data: list[dict[str, object]], lib: str,
          order: list[str], title: str) -> None:
    sub = [row for row in data if row["lib"] == lib]
    x = np.arange(len(order))
    for ind, marker in [("ORT15", "o"), ("ORT16", "^")]:
        by_strategy = {row["strategy"]: row for row in sub if row["ind"] == ind}
        rows = [by_strategy[strategy] for strategy in order]
        ax.plot(x, [row["ts"] for row in rows], marker=marker, ms=8, lw=1.5,
                color=COLOR[ind], label=f"{ind} Ts (C\u2192T/G\u2192A)")
        # transversion baseline (only Untrimmed has data)
        tv = by_strategy["Untrimmed"]["tv"]
        ax.axhline(tv, color=COLOR[ind], ls="--", lw=1, alpha=0.6,
                   label=f"{ind} Tv baseline (untrimmed)")
        # mark Tv point
        ax.plot([0], [tv], marker="s", mfc="white", mec=COLOR[ind],
                ms=8, lw=0)
    # site retention as background bars (right-axis)
    ax2 = ax.twinx()
    width = 0.35
    for k, (ind, off) in enumerate([("ORT15", -width/2), ("ORT16", +width/2)]):
        by_strategy = {row["strategy"]: row for row in sub if row["ind"] == ind}
        retained = np.array([by_strategy[strategy]["retain"] for strategy in order])
        ax2.bar(x + off, retained, width=width,
                color=COLOR[ind], alpha=0.25, edgecolor="none", zorder=0)
        # numeric retention labels above each bar
        for xi, val in zip(x + off, retained):
            if np.isfinite(val):
                ax2.text(xi, val + 1.5, f"{val:.0f}",
                         ha="center", va="bottom", fontsize=7,
                         color=COLOR[ind], alpha=0.9)
    ax2.set_ylim(0, 115)
    ax2.set_ylabel("Sites retained (%)", color="grey", fontsize=9)
    ax2.tick_params(axis="y", colors="grey", labelsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Mean ALT / (REF+ALT) at 1000G sites (%)")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylim(1, 5)
    ax.grid(axis="y", ls=":", alpha=0.4)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = load_data(args.input)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    panel(axes[0], data, "fu", ORDER_FU, "Full-UDG libraries")
    panel(axes[1], data, "nu", ORDER_NU, "Non-UDG libraries")

    # one shared legend
    handles = [
        plt.Line2D([0],[0], marker="o", color=COLOR["ORT15"], lw=1.5, ms=8,
                   label="ORT15  transitions (C\u2192T / G\u2192A)"),
        plt.Line2D([0],[0], marker="^", color=COLOR["ORT16"], lw=1.5, ms=8,
                   label="ORT16  transitions (C\u2192T / G\u2192A)"),
        plt.Line2D([0],[0], marker="s", mfc="white", mec="black",
                   color="black", lw=0, ms=8,
                   label="Transversion baseline (untrimmed)"),
        plt.Line2D([0],[0], color="grey", ls="--", lw=1,
                   label="Tv baseline reference line"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        "Effect of PMD-handling strategies on alternative-allele fraction "
        "at 1000G SNP sites",
        fontsize=13, y=0.99,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
