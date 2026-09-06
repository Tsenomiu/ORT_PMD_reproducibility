#!/usr/bin/env python3
"""Generate the IBD summary and length-versus-count figures.

The script reads the transition-inclusive AADR v62 CROSS64 summary. Figure 4
combines the preselected asymmetric subset (full-UDG held uncorrected,
non-UDG correction varied) with same-method full-UDG/full-UDG and
non-UDG/non-UDG comparisons. Figure S8 is the 16-point length-versus-count
diagnostic across the four library-pair categories.

Primary calls use strict >12 cM (and, for Figure S8 panel b, >20 cM) together
with the >220 SNP/cM density filter. No external gap merge is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


METHODS: Tuple[str, ...] = ("raw", "trim5", "rescale5", "bamrefine5")
METHOD_LABELS: Mapping[str, str] = {
    "raw": "Uncorrected",
    "trim5": "Trim-5",
    "rescale5": "Rescale-5",
    "bamrefine5": "bamRefine-5",
}
METHOD_COLORS: Mapping[str, str] = {
    "raw": "#8F8F8F",
    "trim5": "#C5A15A",
    "rescale5": "#4C78A8",
    "bamrefine5": "#B33A3A",
}
METHOD_MARKERS: Mapping[str, str] = {
    "raw": "o",
    "trim5": "s",
    "rescale5": "D",
    "bamrefine5": "^",
}

PAIR_ORDER: Tuple[str, ...] = (
    "ORT15_fu_x_ORT16_fu",
    "ORT15_fu_x_ORT16_nu",
    "ORT15_nu_x_ORT16_fu",
    "ORT15_nu_x_ORT16_nu",
)
PAIR_LABELS: Mapping[str, str] = {
    "ORT15_fu_x_ORT16_fu": "full-UDG × full-UDG",
    "ORT15_fu_x_ORT16_nu": "full-UDG × non-UDG",
    "ORT15_nu_x_ORT16_fu": "non-UDG × full-UDG",
    "ORT15_nu_x_ORT16_nu": "non-UDG × non-UDG",
}
PAIR_COLORS: Mapping[str, str] = {
    "ORT15_fu_x_ORT16_fu": "#3B6FB6",
    "ORT15_fu_x_ORT16_nu": "#2C8C5A",
    "ORT15_nu_x_ORT16_fu": "#D17A22",
    "ORT15_nu_x_ORT16_nu": "#A83A45",
}

FIGURE4_PANELS = (
    (
        "a",
        "ORT15_fu_x_ORT16_nu",
        "asymmetric",
        "ORT15 full-UDG (uncorrected) × ORT16 non-UDG\n(non-UDG correction varied)",
    ),
    (
        "b",
        "ORT15_nu_x_ORT16_fu",
        "asymmetric",
        "ORT15 non-UDG × ORT16 full-UDG (uncorrected)\n(non-UDG correction varied)",
    ),
    (
        "c",
        "ORT15_fu_x_ORT16_fu",
        "same_treatment",
        "ORT15 full-UDG × ORT16 full-UDG\n(same method on both)",
    ),
    (
        "d",
        "ORT15_nu_x_ORT16_nu",
        "same_treatment",
        "ORT15 non-UDG × ORT16 non-UDG\n(same method on both)",
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def finite_float(row: Mapping[str, str], field: str, context: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid {field!r} in {context}: {row.get(field)!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {field!r} in {context}")
    return value


def nonnegative_integer(row: Mapping[str, str], field: str, context: str) -> int:
    value = finite_float(row, field, context)
    if value < 0 or not value.is_integer():
        raise ValueError(f"Expected non-negative integer {field!r} in {context}")
    return int(value)


def validate_rows(
    rows: Sequence[Mapping[str, str]],
    expected_count: int,
    denominator: float | None = None,
) -> float:
    if len(rows) != expected_count:
        raise ValueError(f"Expected {expected_count} rows, found {len(rows)}")
    observed_denominators = {
        round(finite_float(row, "callable_map_denominator_cM", "summary row"), 6)
        for row in rows
    }
    if len(observed_denominators) != 1:
        raise ValueError(f"Inconsistent callable-map denominators: {observed_denominators}")
    observed = observed_denominators.pop()
    if denominator is not None and abs(observed - denominator) > 1e-5:
        raise ValueError(
            f"Summary denominator {observed:.6f} differs from expected {denominator:.6f}"
        )
    for row in rows:
        context = f"{row.get('iid1')} × {row.get('iid2')}"
        for threshold in (12, 20):
            total = finite_float(row, f"density220_sum_gt{threshold}_cM", context)
            count = nonnegative_integer(row, f"density220_n_gt{threshold}", context)
            recovery = finite_float(row, f"density220_recovery_gt{threshold}_pct", context)
            if total < 0 or count < 0:
                raise ValueError(f"Negative IBD summary in {context}")
            calculated = 100.0 * total / observed
            if abs(calculated - recovery) > 5e-5:
                raise ValueError(
                    f"Recovery mismatch in {context} at >{threshold} cM: "
                    f"reported={recovery:.6f}, calculated={calculated:.6f}"
                )
    return observed


def manifest_denominator(path: Path) -> float:
    with path.open() as handle:
        manifest = json.load(handle)
    source_name = Path(str(manifest.get("v62_source", ""))).name
    if source_name != "v62.0_1240k_public.snp":
        raise ValueError(f"Manifest is not AADR v62.0: {source_name!r}")
    allsnp = manifest.get("marker_sets", {}).get("allsnp", {})
    value = allsnp.get(
        "callable_first_to_last_span_cM",
        manifest.get("canonical_callable_first_to_last_span_cM"),
    )
    denominator = float(value)
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError(f"Invalid manifest denominator: {denominator}")
    return denominator


def index_unique(
    rows: Sequence[Mapping[str, str]], fields: Sequence[str]
) -> Dict[Tuple[str, ...], Mapping[str, str]]:
    indexed: Dict[Tuple[str, ...], Mapping[str, str]] = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        if key in indexed:
            raise ValueError(f"Duplicate design row: {key}")
        indexed[key] = row
    return indexed


def figure4_index(
    asymmetric_rows: Sequence[Mapping[str, str]],
    summary_rows: Sequence[Mapping[str, str]],
) -> Dict[Tuple[str, str], Mapping[str, str]]:
    expected_directions = (
        "ORT15_fu_x_ORT16_nu",
        "ORT15_nu_x_ORT16_fu",
    )
    asymmetric_index = index_unique(
        asymmetric_rows, ("library_direction", "non_udg_treatment")
    )
    expected_keys = {
        (direction, method) for direction in expected_directions for method in METHODS
    }
    if set(asymmetric_index) != expected_keys:
        raise ValueError(
            "Asymmetric design is not the exact two-direction by four-treatment set"
        )

    same_rows = [row for row in summary_rows if row.get("same_treatment") == "yes"]
    same_index = index_unique(
        same_rows, ("library_direction", "same_treatment_method")
    )
    needed_same_keys = {
        (direction, method)
        for direction in ("ORT15_fu_x_ORT16_fu", "ORT15_nu_x_ORT16_nu")
        for method in METHODS
    }
    if not needed_same_keys.issubset(same_index):
        missing = sorted(needed_same_keys.difference(same_index))
        raise ValueError(f"Missing same-treatment Figure 4 rows: {missing}")

    selected: Dict[Tuple[str, str], Mapping[str, str]] = {}
    for panel, direction, design, _ in FIGURE4_PANELS:
        source = asymmetric_index if design == "asymmetric" else same_index
        for method in METHODS:
            selected[(panel, method)] = source[(direction, method)]
    if len(selected) != 16:
        raise ValueError(f"Expected 16 Figure 4 rows, found {len(selected)}")
    return selected


def draw_figure4(
    asymmetric_rows: Sequence[Mapping[str, str]],
    summary_rows: Sequence[Mapping[str, str]],
    output: Path,
) -> Tuple[float, float]:
    indexed = figure4_index(asymmetric_rows, summary_rows)
    all_totals = [
        finite_float(row, "density220_sum_gt12_cM", "Figure 4")
        for row in indexed.values()
    ]
    x_min = math.floor((min(all_totals) - 100.0) / 100.0) * 100.0
    x_max = math.ceil((max(all_totals) + 100.0) / 100.0) * 100.0
    tick_start = math.ceil(x_min / 200.0) * 200.0
    tick_end = math.floor(x_max / 200.0) * 200.0
    ticks = [
        float(value)
        for value in range(int(tick_start), int(tick_end) + 1, 200)
    ]

    figure, axes = plt.subplots(2, 2, figsize=(7.1, 5.55), sharex=True)
    for axis, (panel, direction, _design, panel_title) in zip(
        axes.flat, FIGURE4_PANELS
    ):
        totals = [
            finite_float(indexed[(panel, method)], "density220_sum_gt12_cM", direction)
            for method in METHODS
        ]
        counts = [
            nonnegative_integer(indexed[(panel, method)], "density220_n_gt12", direction)
            for method in METHODS
        ]
        y_positions = list(reversed(range(4)))
        for y_position, method, total, count in zip(
            y_positions, METHODS, totals, counts
        ):
            axis.scatter(
                total,
                y_position,
                s=62,
                color=METHOD_COLORS[method],
                edgecolor="#303030",
                linewidth=0.7,
                zorder=3,
            )
            # Put labels to the left before they can cross into the adjacent panel.
            place_left = total > x_min + 0.55 * (x_max - x_min)
            axis.annotate(
                f"{total:,.0f} cM; $n={count}$",
                xy=(total, y_position),
                xytext=(-7 if place_left else 7, 0),
                textcoords="offset points",
                ha="right" if place_left else "left",
                va="center",
                fontsize=7.2,
            )
        axis.set_yticks(y_positions)
        axis.set_yticklabels([METHOD_LABELS[method] for method in METHODS])
        axis.set_title(panel_title, fontsize=8.2, pad=7, linespacing=1.25)
        axis.set_xlim(x_min, x_max)
        axis.set_xticks(ticks)
        axis.set_ylim(-0.55, 3.55)
        axis.grid(axis="x", color="#D9D9D9", linewidth=0.65)
        axis.set_axisbelow(True)
        axis.text(
            -0.12,
            1.035,
            panel,
            transform=axis.transAxes,
            fontsize=10.5,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
        axis.tick_params(axis="x", labelsize=8.0)
    figure.text(
        0.5,
        0.045,
        "Summed IBD1 length >12 cM (cM)",
        ha="center",
        va="bottom",
        fontsize=8.8,
    )
    figure.text(
        0.5,
        0.012,
        "a-b: full-UDG held uncorrected; c-d: same method on both. Labels give rounded cM and segment count.",
        ha="center",
        va="bottom",
        fontsize=7.0,
        color="#444444",
    )
    figure.subplots_adjust(
        left=0.14, right=0.985, bottom=0.14, top=0.945, wspace=0.31, hspace=0.55
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "Title": "Four-pair IBD1 comparison across PMD correction states",
        "Author": "ORT study",
        "Subject": "Transition-inclusive AADR v62 ancIBD; strict >12 cM and >220 SNP/cM",
        "Creator": "make_figures.py",
    }
    figure.savefig(output, metadata=metadata)
    figure.savefig(output.with_suffix(".png"), dpi=300, metadata={"Software": "Matplotlib"})
    plt.close(figure)
    return x_min, x_max


def draw_length_count_diagnostic(rows: Sequence[Mapping[str, str]], output: Path) -> None:
    selected = [row for row in rows if row.get("same_treatment") == "yes"]
    if len(selected) != 16:
        raise ValueError(f"Expected 16 same-treatment rows, found {len(selected)}")
    indexed = index_unique(selected, ("library_direction", "same_treatment_method"))
    expected_keys = {(pair, method) for pair in PAIR_ORDER for method in METHODS}
    if set(indexed) != expected_keys:
        raise ValueError("Same-treatment design is not the exact four-by-four set")

    figure, axes = plt.subplots(1, 2, figsize=(7.1, 3.55))
    for panel_index, (axis, threshold) in enumerate(zip(axes, (12, 20))):
        for pair in PAIR_ORDER:
            for method in METHODS:
                row = indexed[(pair, method)]
                total = finite_float(row, f"density220_sum_gt{threshold}_cM", pair)
                count = nonnegative_integer(row, f"density220_n_gt{threshold}", pair)
                axis.scatter(
                    total,
                    count,
                    s=47,
                    marker=METHOD_MARKERS[method],
                    facecolor=PAIR_COLORS[pair],
                    edgecolor="#222222",
                    linewidth=0.65,
                    zorder=3,
                )
        axis.set_xlabel(f"Summed IBD1 length >{threshold} cM (cM)", fontsize=8.6)
        axis.set_ylabel(f"Number of segments >{threshold} cM", fontsize=8.6)
        axis.grid(color="#DADADA", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.text(
            -0.14,
            1.03,
            "ab"[panel_index],
            transform=axis.transAxes,
            fontsize=10.5,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    pair_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=PAIR_COLORS[pair],
            markeredgecolor="#222222",
            markersize=6.3,
            label=PAIR_LABELS[pair],
        )
        for pair in PAIR_ORDER
    ]
    method_handles = [
        Line2D(
            [0],
            [0],
            marker=METHOD_MARKERS[method],
            linestyle="none",
            markerfacecolor="#BEBEBE",
            markeredgecolor="#222222",
            markersize=6.3,
            label=METHOD_LABELS[method],
        )
        for method in METHODS
    ]
    figure.legend(
        handles=pair_handles,
        loc="lower left",
        bbox_to_anchor=(0.105, 0.005),
        ncol=2,
        fontsize=7.2,
        title="Library pair (ORT15 × ORT16)",
        title_fontsize=7.4,
        frameon=False,
        handletextpad=0.45,
        columnspacing=1.25,
    )
    figure.legend(
        handles=method_handles,
        loc="lower right",
        bbox_to_anchor=(0.99, 0.005),
        ncol=2,
        fontsize=7.2,
        title="Same method on both libraries",
        title_fontsize=7.4,
        frameon=False,
        handletextpad=0.45,
        columnspacing=1.2,
    )
    figure.subplots_adjust(left=0.105, right=0.985, bottom=0.35, top=0.93, wspace=0.31)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "Title": "Supplementary diagnostic: IBD1 length and segment count",
        "Author": "ORT study",
        "Subject": "AADR v62 all-SNP ancIBD with >220 SNP/cM density filter",
        "Creator": "make_figures.py",
    }
    figure.savefig(output, metadata=metadata)
    figure.savefig(output.with_suffix(".png"), dpi=300, metadata={"Software": "Matplotlib"})
    plt.close(figure)


def write_validation(
    output_dir: Path,
    asymmetric_path: Path,
    summary_path: Path,
    manifest_path: Path,
    denominator: float,
    asymmetric_rows: Sequence[Mapping[str, str]],
    summary_rows: Sequence[Mapping[str, str]],
    figure4_xlim: Tuple[float, float],
) -> None:
    payload = {
        "status": "PASS",
        "analysis": "AADR v62 all-SNP corrected ancIBD",
        "filters": {
            "length": "strict >12 cM; Figure S8 also strict >20 cM",
            "density": "strict >220 SNP/cM",
            "external_gap_merge": False,
        },
        "callable_map_denominator_cM": denominator,
        "inputs": {
            str(asymmetric_path): sha256(asymmetric_path),
            str(summary_path): sha256(summary_path),
            str(manifest_path): sha256(manifest_path),
        },
        "design": {
            "figure4_panels": 4,
            "figure4_rows": len(figure4_index(asymmetric_rows, summary_rows)),
            "figure4_shared_x_limits_cM": list(figure4_xlim),
            "figure4_panel_definitions": [
                {
                    "panel": panel,
                    "library_direction": direction,
                    "design": design,
                }
                for panel, direction, design, _title in FIGURE4_PANELS
            ],
            "figureS8_same_treatment_rows": sum(
                row.get("same_treatment") == "yes" for row in summary_rows
            ),
        },
    }
    with (output_dir / "FIG4_S8_VALIDATION.json").open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_display_values(
    output_dir: Path,
    asymmetric_rows: Sequence[Mapping[str, str]],
    summary_rows: Sequence[Mapping[str, str]],
) -> None:
    selected_index = figure4_index(asymmetric_rows, summary_rows)
    with (output_dir / "FIG4_DISPLAY_VALUES.tsv").open("w", newline="") as handle:
        fields = (
            "panel",
            "design",
            "library_direction",
            "method",
            "sum_gt12_cM_exact",
            "sum_gt12_cM_display",
            "n_gt12",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for panel, direction, design, _title in FIGURE4_PANELS:
            for method in METHODS:
                row = selected_index[(panel, method)]
                total = finite_float(row, "density220_sum_gt12_cM", direction)
                writer.writerow(
                    {
                        "panel": panel,
                        "design": design,
                        "library_direction": direction,
                        "method": method,
                        "sum_gt12_cM_exact": f"{total:.6f}",
                        "sum_gt12_cM_display": f"{total:.0f}",
                        "n_gt12": nonnegative_integer(
                            row, "density220_n_gt12", direction
                        ),
                    }
                )

    selected = [row for row in summary_rows if row.get("same_treatment") == "yes"]
    same_index = index_unique(
        selected, ("library_direction", "same_treatment_method")
    )
    with (output_dir / "FIGS8_DISPLAY_VALUES.tsv").open("w", newline="") as handle:
        fields = (
            "library_direction",
            "same_treatment_method",
            "sum_gt12_cM_exact",
            "n_gt12",
            "sum_gt20_cM_exact",
            "n_gt20",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for pair in PAIR_ORDER:
            for method in METHODS:
                row = same_index[(pair, method)]
                writer.writerow(
                    {
                        "library_direction": pair,
                        "same_treatment_method": method,
                        "sum_gt12_cM_exact": f"{finite_float(row, 'density220_sum_gt12_cM', pair):.6f}",
                        "n_gt12": nonnegative_integer(
                            row, "density220_n_gt12", pair
                        ),
                        "sum_gt20_cM_exact": f"{finite_float(row, 'density220_sum_gt20_cM', pair):.6f}",
                        "n_gt20": nonnegative_integer(
                            row, "density220_n_gt20", pair
                        ),
                    }
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asymmetric", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--resource-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    denominator = manifest_denominator(args.resource_manifest)
    _, asymmetric_rows = read_tsv(args.asymmetric)
    _, summary_rows = read_tsv(args.summary)
    validate_rows(asymmetric_rows, 8, denominator)
    validate_rows(summary_rows, 64, denominator)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure4 = args.output_dir / "figure_04_ibd_summary.pdf"
    figure_s8 = args.output_dir / "supplementary_08_ibd_length_count.pdf"
    figure4_xlim = draw_figure4(asymmetric_rows, summary_rows, figure4)
    draw_length_count_diagnostic(summary_rows, figure_s8)
    write_display_values(args.output_dir, asymmetric_rows, summary_rows)
    write_validation(
        args.output_dir,
        args.asymmetric,
        args.summary,
        args.resource_manifest,
        denominator,
        asymmetric_rows,
        summary_rows,
        figure4_xlim,
    )
    print(f"figure4={figure4}")
    print(f"figureS8={figure_s8}")
    print(f"callable_map_denominator_cM={denominator:.6f}")
    print("design_validation=PASS (4 Figure 4 panels/16 rows; 16 Figure S8 rows)")


if __name__ == "__main__":
    main()
