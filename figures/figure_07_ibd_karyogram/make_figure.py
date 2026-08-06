#!/usr/bin/env python3
"""Plot deterministic 4x2 ancIBD karyograms from the corrected v62 run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


METHODS = ("raw", "trim5", "rescale5", "bamrefine5")
METHOD_LABELS = {
    "raw": "Raw",
    "trim5": "Trim-5",
    "rescale5": "Rescale-5",
    "bamrefine5": "bamRefine-5",
}
DIRECTIONS = ("fu_x_nu", "nu_x_fu")
DIRECTION_LABELS = {
    "fu_x_nu": "ORT15 full-UDG × ORT16 non-UDG",
    "nu_x_fu": "ORT15 non-UDG × ORT16 full-UDG",
}


def marker_set_denominator(manifest: dict, marker_set: str) -> float:
    metadata = manifest.get("marker_sets", {}).get(marker_set)
    if metadata:
        return float(metadata["callable_first_to_last_span_cM"])
    if marker_set == "transversion":
        return float(manifest["transversion_callable_first_to_last_span_cM"])
    return float(
        manifest.get(
            "canonical_callable_first_to_last_span_cM",
            manifest["callable_first_to_last_span_cM"],
        )
    )


def chromosome_bounds(manifest: dict, marker_set: str, chromosome: int):
    record = manifest["chromosomes"][str(chromosome)]
    if marker_set == "transversion":
        return (
            float(record["transversion_canonical_first_M"]),
            float(record["transversion_canonical_last_M"]),
        )
    first = record["canonical_first_M"] if "canonical_first_M" in record else record["map_first_M"]
    last = record["canonical_last_M"] if "canonical_last_M" in record else record["map_last_M"]
    return (float(first), float(last))


def check_summary_marker_set(summaries, marker_set: str) -> None:
    observed = {row.get("marker_set", "allsnp") for row in summaries}
    if observed != {marker_set}:
        raise ValueError(
            f"Summary marker set {sorted(observed)} does not match requested {marker_set}"
        )


def check_base_marker_set(base: Path, marker_set: str) -> None:
    path = base / "MARKER_SET.txt"
    if not path.is_file():
        if marker_set == "transversion":
            raise ValueError("Transversion output is missing MARKER_SET.txt")
        return
    metadata = dict(
        raw.split("=", 1)
        for raw in path.read_text().splitlines()
        if "=" in raw
    )
    if metadata.get("marker_set") != marker_set:
        raise ValueError(
            f"Output tree marker set {metadata.get('marker_set')!r} does not "
            f"match requested {marker_set!r}"
        )


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def plot_mode(
    base: Path, manifest: dict, output: Path, mode: str, marker_set: str = "allsnp"
) -> None:
    denominator = marker_set_denominator(manifest, marker_set)
    panel_label = "all-SNP" if marker_set == "allsnp" else "transversion-only"
    if mode == "density220":
        segment_path = base / "FIG8_SEGMENTS_DENSITY220.tsv"
        sum_prefix = "density220"
        title = (
            f"Corrected {panel_label} ancIBD: native IBD1 segments passing "
            "the density filter"
        )
        subtitle = (
            "Exact AADR v62/AF/GLIMPSE intersection; segments >8 cM and "
            ">220 SNP/cM; "
            "no external 6-cM merge"
        )
        color = "#4C78A8"
    elif mode == "native":
        segment_path = base / "FIG8_SEGMENTS_NATIVE.tsv"
        sum_prefix = "native"
        title = f"Corrected {panel_label} ancIBD: native IBD1 segment calls"
        subtitle = (
            "Exact AADR v62/AF/GLIMPSE intersection; segments >8 cM; "
            "no SNP-density filter and "
            "no external 6-cM merge"
        )
        color = "#B279A2"
    else:
        raise ValueError(mode)

    segments = read_tsv(segment_path)
    summaries = read_tsv(base / "FIG8_PAIR_SUMMARY.tsv")
    check_summary_marker_set(summaries, marker_set)
    summary_index = {(row["method"], row["direction"]): row for row in summaries}

    starts = {
        chromosome: chromosome_bounds(manifest, marker_set, chromosome)[0] * 100
        for chromosome in range(1, 23)
    }
    spans = {
        chromosome: (
            chromosome_bounds(manifest, marker_set, chromosome)[1]
            - chromosome_bounds(manifest, marker_set, chromosome)[0]
        )
        * 100
        for chromosome in range(1, 23)
    }

    figure, axes = plt.subplots(4, 2, figsize=(8.2, 9.7), sharex=True, sharey=True)
    background = "#D9D9D9"

    for row_index, method in enumerate(METHODS):
        for column_index, direction in enumerate(DIRECTIONS):
            axis = axes[row_index, column_index]
            selected = [
                row
                for row in segments
                if row["method"] == method and row["direction"] == direction
            ]

            for chromosome in range(1, 23):
                y = 22 - chromosome
                axis.add_patch(
                    Rectangle(
                        (0, y - 0.31),
                        spans[chromosome],
                        0.62,
                        facecolor=background,
                        edgecolor="none",
                        zorder=1,
                    )
                )

            for segment in selected:
                chromosome = int(segment["chromosome"])
                start = float(segment["start_cM"]) - starts[chromosome]
                end = float(segment["end_cM"]) - starts[chromosome]
                start = max(0.0, min(spans[chromosome], start))
                end = max(0.0, min(spans[chromosome], end))
                if end <= start:
                    continue
                y = 22 - chromosome
                axis.add_patch(
                    Rectangle(
                        (start, y - 0.31),
                        end - start,
                        0.62,
                        facecolor=color,
                        edgecolor="none",
                        zorder=2,
                    )
                )

            summary = summary_index[(method, direction)]
            summary_denominator = float(summary["callable_map_denominator_cM"])
            if abs(summary_denominator - denominator) > 0.001:
                raise ValueError(
                    "Figure summary and resource manifest use different map denominators: "
                    f"{summary_denominator} vs {denominator}"
                )
            total = float(summary[f"{sum_prefix}_sum_gt8_cM"])
            count = int(summary[f"{sum_prefix}_n_gt8"])
            axis.text(
                0.98,
                0.035,
                f"Sum IBD1 >8 cM: {total:.0f} cM · n={count}",
                transform=axis.transAxes,
                ha="right",
                va="bottom",
                fontsize=7.2,
                fontweight="bold",
                color=color,
            )

            axis.set_xlim(0, 300)
            axis.set_ylim(-0.7, 22.2)
            axis.set_xticks((0, 50, 100, 150, 200, 250, 300))
            if column_index == 0:
                axis.set_ylabel(METHOD_LABELS[method], fontsize=8)
                ticks = [22 - chromosome for chromosome in (1, 5, 10, 15, 20)]
                axis.set_yticks(ticks)
                axis.set_yticklabels((1, 5, 10, 15, 20), fontsize=6)
            if row_index == 0:
                axis.set_title(DIRECTION_LABELS[direction], fontsize=8)
            if row_index == len(METHODS) - 1:
                axis.set_xlabel("Position from first v62 marker (cM)", fontsize=7)
            axis.tick_params(axis="x", labelsize=6)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)

    figure.suptitle(title, fontsize=10, y=0.992)
    figure.text(0.5, 0.957, subtitle, ha="center", va="top", fontsize=7, style="italic")
    figure.text(0.015, 0.5, "Chromosome", va="center", rotation=90, fontsize=7)
    figure.text(
        0.5,
        0.012,
        (
            "Recovery denominator: first-to-last autosomal span of the exact "
            f"v62/AF/GLIMPSE intersection ({denominator:.2f} cM)."
        ),
        ha="center",
        fontsize=6.5,
        color="#555555",
    )
    figure.tight_layout(rect=(0.035, 0.03, 1, 0.94))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    figure.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_ancibd_style(
    base: Path,
    manifest: dict,
    output: Path,
    threshold_cm: float = 12.0,
    marker_set: str = "allsnp",
    segment_mode: str = "density220",
) -> None:
    """Plot IBD1 calls using ancIBD's two-row chromosome geometry.

    The installed ancIBD 0.7 helper represents each chromosome as a rounded
    vertical bar, with chromosomes 1--11 above chromosomes 12--22.  This
    function retains that geometry while arranging the eight predefined
    treatment comparisons in one deterministic 4-by-2 panel.
    """

    if threshold_cm not in (8.0, 12.0, 16.0, 20.0):
        raise ValueError("Threshold must match a collected strict ancIBD bin")
    if segment_mode not in ("native", "density220"):
        raise ValueError(segment_mode)
    if marker_set == "transversion" and segment_mode != "native":
        raise ValueError(
            "Transversion inputs require native calls rather than the 220-SNP/cM filter"
        )

    denominator = marker_set_denominator(manifest, marker_set)
    segment_filename = (
        "FIG8_SEGMENTS_NATIVE.tsv"
        if segment_mode == "native"
        else "FIG8_SEGMENTS_DENSITY220.tsv"
    )
    segments = read_tsv(base / segment_filename)
    summaries = read_tsv(base / "FIG8_PAIR_SUMMARY.tsv")
    check_summary_marker_set(summaries, marker_set)
    summary_index = {(row["method"], row["direction"]): row for row in summaries}
    starts_m = {
        chromosome: chromosome_bounds(manifest, marker_set, chromosome)[0]
        for chromosome in range(1, 23)
    }
    spans_m = {
        chromosome: chromosome_bounds(manifest, marker_set, chromosome)[1]
        - starts_m[chromosome]
        for chromosome in range(1, 23)
    }

    # Use ancIBD's native visual grammar: separate chromosome rows with a
    # Morgan y-scale, chromosomes 1--11 above 12--22, rounded gray backbones,
    # and flat-ended colored IBD intervals.  A slightly taller canvas and
    # lighter, outline-free backbones keep small uncalled gaps visible without
    # turning the eight-panel comparison into a heavy barcode.
    figure = plt.figure(figsize=(7.1, 7.6))
    outer_grid = gridspec.GridSpec(
        4,
        2,
        figure=figure,
        left=0.105,
        right=0.995,
        bottom=0.12,
        top=0.94,
        hspace=1.02,
        wspace=0.17,
    )
    chromosome_color = "#CECECE"
    ibd_color = "#443A8B"
    threshold_key = int(threshold_cm)

    for row_index, method in enumerate(METHODS):
        for column_index, direction in enumerate(DIRECTIONS):
            inner_grid = gridspec.GridSpecFromSubplotSpec(
                2,
                1,
                subplot_spec=outer_grid[row_index, column_index],
                height_ratios=(3, 2),
                hspace=0.62,
            )
            upper_axis = figure.add_subplot(inner_grid[0])
            lower_axis = figure.add_subplot(inner_grid[1])
            selected = [
                row
                for row in segments
                if row["method"] == method
                and row["direction"] == direction
                and float(row["length_cM"]) > threshold_cm
            ]

            for chromosome in range(1, 23):
                x_position = chromosome if chromosome <= 11 else chromosome - 11
                axis = upper_axis if chromosome <= 11 else lower_axis
                axis.plot(
                    [x_position, x_position],
                    [-0.03, spans_m[chromosome] + 0.03],
                    linewidth=6.0,
                    color=chromosome_color,
                    solid_capstyle="round",
                    zorder=0,
                    clip_on=False,
                )

            for segment in selected:
                chromosome = int(segment["chromosome"])
                x_position = chromosome if chromosome <= 11 else chromosome - 11
                axis = upper_axis if chromosome <= 11 else lower_axis
                start_m = float(segment["start_cM"]) / 100.0 - starts_m[chromosome]
                end_m = float(segment["end_cM"]) / 100.0 - starts_m[chromosome]
                start_m = max(0.0, min(spans_m[chromosome], start_m))
                end_m = max(0.0, min(spans_m[chromosome], end_m))
                if end_m <= start_m:
                    continue
                axis.plot(
                    [x_position, x_position],
                    [start_m, end_m],
                    linewidth=5.0,
                    color=ibd_color,
                    solid_capstyle="butt",
                    zorder=1,
                )

            summary = summary_index[(method, direction)]
            total = float(summary[f"{segment_mode}_sum_gt{threshold_key}_cM"])
            lower_axis.text(
                0.50,
                -0.58,
                f"{total:,.0f} cM",
                transform=lower_axis.transAxes,
                ha="center",
                va="top",
                fontsize=8.1,
                fontweight="semibold",
                color="#2F2F2F",
                clip_on=False,
            )

            upper_axis.set_xlim(0.3, 11.5)
            lower_axis.set_xlim(0.3, 11.5)
            upper_axis.set_ylim(-0.22, 3.20)
            lower_axis.set_ylim(-0.22, 2.05)
            upper_axis.set_xticks(range(1, 12))
            upper_axis.set_xticklabels(range(1, 12), fontsize=8.0)
            lower_axis.set_xticks(range(1, 12))
            lower_axis.set_xticklabels(range(12, 23), fontsize=8.0)

            for axis in (upper_axis, lower_axis):
                axis.spines["right"].set_visible(False)
                axis.spines["top"].set_visible(False)
                axis.spines["bottom"].set_visible(False)
                axis.tick_params(axis="x", length=0, pad=2.5)
                axis.tick_params(axis="y", labelsize=8.0, length=2.2, pad=1)
                if column_index == 1:
                    axis.spines["left"].set_visible(False)
                    axis.tick_params(axis="y", left=False, labelleft=False)

            upper_axis.set_yticks((0, 1, 2, 3))
            lower_axis.set_yticks((0, 1))
            if column_index == 0:
                lower_axis.text(
                    -0.105,
                    1.34,
                    METHOD_LABELS[method],
                    transform=lower_axis.transAxes,
                    ha="center",
                    va="center",
                    rotation=90,
                    fontsize=8.5,
                    fontweight="medium",
                    clip_on=False,
                )
            if row_index == 0:
                upper_axis.set_title(
                    DIRECTION_LABELS[direction],
                    fontsize=9.4,
                    fontweight="medium",
                    pad=5,
                )

    figure.text(
        0.020,
        0.49,
        "Genetic map span (M)",
        rotation=90,
        ha="center",
        va="center",
        fontsize=8.0,
    )
    figure.text(0.055, 0.962, "a", fontsize=10.5, fontweight="bold")
    figure.text(0.535, 0.962, "b", fontsize=10.5, fontweight="bold")
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=chromosome_color,
            linewidth=6.0,
            solid_capstyle="round",
        ),
        Line2D([0], [0], color=ibd_color, linewidth=5.0, solid_capstyle="butt"),
    ]
    figure.legend(
        legend_handles,
        (
            "Callable autosomal span",
            f"IBD1 segment >{threshold_key} cM; >220 SNP/cM",
        ),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.008),
        ncol=2,
        frameon=False,
        fontsize=8.2,
        handlelength=2.1,
        columnspacing=1.8,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight", facecolor="white")
    figure.savefig(
        output.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white"
    )
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--resource-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.resource_manifest.open() as handle:
        manifest = json.load(handle)
    check_base_marker_set(args.base, "allsnp")
    plot_ancibd_style(
        args.base,
        manifest,
        args.output,
        threshold_cm=12.0,
        marker_set="allsnp",
        segment_mode="density220",
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
