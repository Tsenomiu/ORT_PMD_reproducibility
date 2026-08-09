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

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


METHODS = ("raw", "trim5", "rescale5", "bamrefine5")
METHOD_LABELS = {
    "raw": "Uncorrected",
    "trim5": "Trim-5",
    "rescale5": "Rescale-5",
    "bamrefine5": "bamRefine-5",
}
DIRECTIONS = ("fu_x_nu", "nu_x_fu")
DIRECTION_LABELS = {
    "fu_x_nu": "ORT15 full-UDG × ORT16 non-UDG",
    "nu_x_fu": "ORT15 non-UDG × ORT16 full-UDG",
}
CHROMOSOME_X_SPACING = 0.76


def chromosome_x_position(chromosome: int) -> float:
    """Return a centred, slightly compact x-position within either row."""

    row_index = chromosome if chromosome <= 11 else chromosome - 11
    return 6.0 + (row_index - 6.0) * CHROMOSOME_X_SPACING


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


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def plot_ancibd_style(
    base: Path,
    manifest: dict,
    output: Path,
    threshold_cm: float = 12.0,
    marker_set: str = "allsnp",
    segment_mode: str = "density220",
) -> None:
    """Plot the corrected calls using ancIBD's two-row chromosome geometry.

    The installed ancIBD 0.7 helper represents each chromosome as a rounded
    vertical bar, with chromosomes 1--11 above chromosomes 12--22.  This
    function retains that geometry while arranging the eight predefined
    treatment comparisons in one deterministic 4-by-2 manuscript panel.
    """

    if threshold_cm not in (8.0, 12.0, 16.0, 20.0):
        raise ValueError("Threshold must match a collected strict ancIBD bin")
    if marker_set != "allsnp" or segment_mode != "density220":
        raise ValueError("The public Figure 5 generator uses the all-SNP density-filtered input")

    denominator = marker_set_denominator(manifest, marker_set)
    segment_filename = "KARYOGRAM_SEGMENTS_DENSITY220.tsv"
    segments = read_tsv(base / segment_filename)
    summaries = read_tsv(base / "KARYOGRAM_PAIR_SUMMARY.tsv")
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
    # and flat-ended colored IBD intervals.  The compact canvas and reduced
    # inter-row spacing let the manuscript use nearly the full text width,
    # making chromosomes and short retained intervals easier to read without
    # changing their genetic-map geometry.
    figure = plt.figure(figsize=(7.2, 6.55))
    outer_grid = gridspec.GridSpec(
        4,
        2,
        figure=figure,
        left=0.105,
        right=0.995,
        bottom=0.105,
        top=0.94,
        hspace=0.55,
        wspace=0.17,
    )
    chromosome_color = "#D5D5D5"
    ibd_color = "#443A8B"
    threshold_key = int(threshold_cm)

    for row_index, method in enumerate(METHODS):
        for column_index, direction in enumerate(DIRECTIONS):
            inner_grid = gridspec.GridSpecFromSubplotSpec(
                2,
                1,
                subplot_spec=outer_grid[row_index, column_index],
                height_ratios=(3, 2),
                hspace=0.34,
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
                x_position = chromosome_x_position(chromosome)
                axis = upper_axis if chromosome <= 11 else lower_axis
                axis.plot(
                    [x_position, x_position],
                    [-0.03, spans_m[chromosome] + 0.03],
                    linewidth=8.0,
                    color=chromosome_color,
                    solid_capstyle="round",
                    zorder=0,
                    clip_on=False,
                )

            for segment in selected:
                chromosome = int(segment["chromosome"])
                x_position = chromosome_x_position(chromosome)
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
                    linewidth=6.6,
                    color=ibd_color,
                    solid_capstyle="butt",
                    zorder=1,
                )

            summary = summary_index[(method, direction)]
            total = float(summary[f"{segment_mode}_sum_gt{threshold_key}_cM"])
            lower_axis.text(
                0.50,
                -0.52,
                f"{total:,.0f} cM",
                transform=lower_axis.transAxes,
                ha="center",
                va="top",
                fontsize=8.3,
                fontweight="semibold",
                color="#2F2F2F",
                clip_on=False,
            )

            upper_axis.set_xlim(0.3, 11.5)
            lower_axis.set_xlim(0.3, 11.5)
            upper_axis.set_ylim(-0.22, 3.20)
            lower_axis.set_ylim(-0.22, 2.05)
            upper_axis.set_xticks(
                [chromosome_x_position(chromosome) for chromosome in range(1, 12)]
            )
            upper_axis.set_xticklabels(range(1, 12), fontsize=8.3)
            lower_axis.set_xticks(
                [chromosome_x_position(chromosome) for chromosome in range(12, 23)]
            )
            lower_axis.set_xticklabels(range(12, 23), fontsize=8.3)

            for axis in (upper_axis, lower_axis):
                axis.spines["right"].set_visible(False)
                axis.spines["top"].set_visible(False)
                axis.spines["bottom"].set_visible(False)
                axis.tick_params(axis="x", length=0, pad=2.5)
                axis.tick_params(axis="y", labelsize=8.3, length=2.2, pad=1)
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
                    fontsize=8.8,
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
        fontsize=8.3,
    )
    figure.text(0.055, 0.962, "a", fontsize=10.5, fontweight="bold")
    figure.text(0.535, 0.962, "b", fontsize=10.5, fontweight="bold")
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=chromosome_color,
            linewidth=8.0,
            solid_capstyle="round",
        ),
        Line2D([0], [0], color=ibd_color, linewidth=6.6, solid_capstyle="butt"),
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
        fontsize=8.4,
        handlelength=2.1,
        columnspacing=1.8,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output,
        bbox_inches="tight",
        facecolor="white",
        metadata={"CreationDate": None},
    )
    figure.savefig(
        output.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "ORT PMD deterministic figure generator"},
    )
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--resource-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    with args.resource_manifest.open() as handle:
        manifest = json.load(handle)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    marker_set = "allsnp"
    output_name = "figure_05_ibd_karyogram.pdf"
    plot_ancibd_style(
        args.base,
        manifest,
        args.output_dir / output_name,
        threshold_cm=12.0,
        marker_set=marker_set,
        segment_mode="density220",
    )
    caption = (
        "Native ancIBD IBD1 calls from the exact AADR v62/AF/common-GLIMPSE "
        "all-SNP intersection. Segments are shown only when they are strictly "
        "longer than 12 cM and exceed 220 SNP/cM. No external gap merging was "
        "applied."
    )

    marker_metadata = manifest.get("marker_sets", {}).get(marker_set, {})
    marker_count = marker_metadata.get("marker_count")
    if marker_count is None:
        marker_count = manifest.get("totals", {}).get("exact")
    mean_density = marker_metadata.get("mean_marker_density_per_cM")
    if mean_density is None and marker_count is not None:
        mean_density = float(marker_count) / marker_set_denominator(
            manifest, marker_set
        )

    resource_manifest_sha = sha256(args.resource_manifest)
    resource_set_sha = marker_metadata.get(
        "resource_set_sha256", manifest.get("resource_set_sha256")
    )
    if resource_set_sha is None:
        resource_set_sha = resource_manifest_sha

    plot_metadata = {
        "marker_set": marker_set,
        "primary_summary": "density220",
        "length_threshold_semantics": "strict_greater_than",
        "external_gap_merge": False,
        "callable_first_to_last_span_cM": marker_set_denominator(
            manifest, marker_set
        ),
        "marker_count": marker_count,
        "mean_marker_density_per_cM": mean_density,
        "resource_manifest_sha256": resource_manifest_sha,
        "resource_set_sha256": resource_set_sha,
        "caption_template": caption,
        "outputs": [output_name],
    }
    if manifest.get("source_resource_manifest_sha256") is not None:
        plot_metadata["source_resource_manifest_sha256"] = manifest[
            "source_resource_manifest_sha256"
        ]
    (args.output_dir / "figure_05_plot_metadata.json").write_text(
        json.dumps(plot_metadata, indent=2, sort_keys=True) + "\n"
    )
    print(f"marker_set={marker_set}")
    print(f"output_dir={args.output_dir}")


if __name__ == "__main__":
    main()
