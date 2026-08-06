#!/usr/bin/env python3
"""Plot validated 8-by-8 cross-treatment ancIBD summary heatmaps.

This is a deterministic custom visualization of values already collected in
``CROSS64_PAIR_SUMMARY.tsv``.  It is not an ancIBD-native plotting helper.
The primary view uses strict >12 cM and >220 SNP/cM filters; a second view
shows the strict >12 cM native calls without a SNP-density filter.
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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle


matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    }
)


METHODS: Tuple[str, ...] = ("raw", "trim5", "rescale5", "bamrefine5")
METHOD_LABELS: Mapping[str, str] = {
    "raw": "Raw",
    "trim5": "Trim-5",
    "rescale5": "Rescale-5",
    "bamrefine5": "bamRefine-5",
}
LIBRARIES: Tuple[str, ...] = ("fu", "nu")
LIBRARY_LABELS: Mapping[str, str] = {
    "fu": "full-UDG",
    "nu": "non-UDG",
}

# Deliberately truncated, explicitly labelled scales make the modest
# cross-treatment contrasts visible. Exact totals and segment counts remain
# printed in every cell, so the color mapping cannot substitute for the data.
COLOR_LIMITS_CM: Mapping[str, Tuple[float, float]] = {
    "density220": (2700.0, 3300.0),
    "native": (2900.0, 3500.0),
}


def sample_order(individual: str) -> List[str]:
    return [
        f"{individual}_{library}_{method}"
        for library in LIBRARIES
        for method in METHODS
    ]


ORT15_ORDER = sample_order("ORT15")
ORT16_ORDER = sample_order("ORT16")
EXPECTED_PAIRS = [
    (iid1, iid2) for iid1 in ORT15_ORDER for iid2 in ORT16_ORDER
]


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
        raise ValueError(
            f"Invalid numeric field {field!r} in {context}: {row.get(field)!r}"
        ) from error
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {field!r} in {context}: {value}")
    return value


def integral_count(row: Mapping[str, str], field: str, context: str) -> int:
    value = finite_float(row, field, context)
    if value < 0 or not value.is_integer():
        raise ValueError(f"Expected a non-negative integer {field!r} in {context}")
    return int(value)


def manifest_denominator(manifest: Mapping) -> float:
    v62_source = str(manifest.get("v62_source", ""))
    if Path(v62_source).name != "v62.0_1240k_public.snp":
        raise ValueError(
            "Resource manifest is not bound to AADR v62.0: "
            f"v62_source={v62_source!r}"
        )

    marker_sets = manifest.get("marker_sets", {})
    allsnp = marker_sets.get("allsnp", {}) if isinstance(marker_sets, dict) else {}
    if "callable_first_to_last_span_cM" in allsnp:
        denominator = float(allsnp["callable_first_to_last_span_cM"])
    elif "canonical_callable_first_to_last_span_cM" in manifest:
        denominator = float(manifest["canonical_callable_first_to_last_span_cM"])
    else:
        raise ValueError(
            "Resource manifest lacks the exact all-SNP callable-map denominator"
        )
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError(f"Invalid callable-map denominator: {denominator}")
    return denominator


def split_iid(iid: str, expected_individual: str) -> Tuple[str, str]:
    prefix = expected_individual + "_"
    if not iid.startswith(prefix):
        raise ValueError(f"Expected {expected_individual} IID, found {iid!r}")
    fields = iid[len(prefix) :].split("_", 1)
    if len(fields) != 2:
        raise ValueError(f"Malformed IID: {iid!r}")
    library, method = fields
    if library not in LIBRARIES or method not in METHODS:
        raise ValueError(f"Unexpected library/treatment in IID: {iid!r}")
    return library, method


def validate_metadata(row: Mapping[str, str], row_number: int) -> None:
    iid1 = row["iid1"]
    iid2 = row["iid2"]
    library1, treatment1 = split_iid(iid1, "ORT15")
    library2, treatment2 = split_iid(iid2, "ORT16")
    expected = {
        "individual1": "ORT15",
        "individual2": "ORT16",
        "library1": library1,
        "library2": library2,
        "treatment1": treatment1,
        "treatment2": treatment2,
        "library_direction": f"ORT15_{library1}_x_ORT16_{library2}",
        "same_treatment": "yes" if treatment1 == treatment2 else "no",
        "mixed_library": "yes" if library1 != library2 else "no",
        "same_treatment_method": treatment1 if treatment1 == treatment2 else "",
        "non_udg_treatment": (
            treatment1
            if library1 == "nu" and library2 == "fu"
            else treatment2 if library1 == "fu" and library2 == "nu" else ""
        ),
    }
    for field, expected_value in expected.items():
        observed = row.get(field)
        if observed != expected_value:
            raise ValueError(
                f"Row {row_number} has inconsistent {field}: "
                f"observed={observed!r}, expected={expected_value!r}"
            )


def validate_and_build_matrices(
    summary_path: Path,
    manifest_path: Path,
) -> Tuple[float, Dict[str, List[List[float]]], Dict[str, List[List[int]]]]:
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    denominator = manifest_denominator(manifest)

    fields, rows = read_tsv(summary_path)
    required = {
        "iid1",
        "iid2",
        "individual1",
        "individual2",
        "library1",
        "library2",
        "treatment1",
        "treatment2",
        "library_direction",
        "same_treatment",
        "mixed_library",
        "same_treatment_method",
        "non_udg_treatment",
        "callable_map_denominator_cM",
        "density220_sum_gt12_cM",
        "density220_n_gt12",
        "density220_recovery_gt12_pct",
        "native_sum_gt12_cM",
        "native_n_gt12",
        "native_recovery_gt12_pct",
    }
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"CROSS64 summary lacks required columns: {missing}")
    if len(rows) != 64:
        raise ValueError(f"Expected exactly 64 summary rows, found {len(rows)}")

    observed_pairs = [(row["iid1"], row["iid2"]) for row in rows]
    if observed_pairs != EXPECTED_PAIRS:
        first_mismatch = next(
            (
                index
                for index, (observed, expected) in enumerate(
                    zip(observed_pairs, EXPECTED_PAIRS), start=1
                )
                if observed != expected
            ),
            None,
        )
        raise ValueError(
            "CROSS64 summary is not the exact ordered 8-by-8 ORT15-by-ORT16 "
            f"design; first mismatch row={first_mismatch}"
        )
    if len(set(observed_pairs)) != 64:
        raise ValueError("CROSS64 summary contains duplicate pairs")

    matrices: Dict[str, List[List[float]]] = {
        "density220": [[0.0 for _ in range(8)] for _ in range(8)],
        "native": [[0.0 for _ in range(8)] for _ in range(8)],
    }
    counts: Dict[str, List[List[int]]] = {
        "density220": [[0 for _ in range(8)] for _ in range(8)],
        "native": [[0 for _ in range(8)] for _ in range(8)],
    }

    for index, row in enumerate(rows):
        row_number = index + 2
        context = f"{summary_path.name}:{row_number}"
        validate_metadata(row, row_number)
        observed_denominator = finite_float(
            row, "callable_map_denominator_cM", context
        )
        if abs(observed_denominator - denominator) > 1e-5:
            raise ValueError(
                f"{context} denominator {observed_denominator:.9f} does not "
                f"match resource manifest {denominator:.9f}"
            )

        matrix_row = index // 8
        matrix_column = index % 8
        for mode in ("density220", "native"):
            total = finite_float(row, f"{mode}_sum_gt12_cM", context)
            count = integral_count(row, f"{mode}_n_gt12", context)
            reported = finite_float(row, f"{mode}_recovery_gt12_pct", context)
            calculated = 100.0 * total / denominator
            if total < 0:
                raise ValueError(f"Negative IBD1 total in {context}")
            if calculated < -1e-8 or calculated > 100.0 + 1e-6:
                raise ValueError(
                    f"{mode} called map fraction outside 0-100% in {context}: "
                    f"{calculated:.8f}%"
                )
            if abs(reported - calculated) > 5e-5:
                raise ValueError(
                    f"{mode} reported recovery disagrees with sum/denominator in "
                    f"{context}: reported={reported:.6f}, calculated={calculated:.6f}"
                )
            matrices[mode][matrix_row][matrix_column] = calculated
            counts[mode][matrix_row][matrix_column] = count

    return denominator, matrices, counts


def annotation_color(rgba: Sequence[float]) -> str:
    """Choose black or white text using WCAG relative contrast.

    Matplotlib returns sRGB channel values.  WCAG relative luminance requires
    linearising those channels before applying the luminance weights; using
    the weights directly on sRGB can leave mid-tone viridis cells with poor
    label contrast.
    """

    def linearise(channel: float) -> float:
        return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4

    red, green, blue = (linearise(float(channel)) for channel in rgba[:3])
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    contrast_black = (luminance + 0.05) / 0.05
    contrast_white = 1.05 / (luminance + 0.05)
    return "#000000" if contrast_black >= contrast_white else "#FFFFFF"


def draw_heatmap(
    values: Sequence[Sequence[float]],
    counts: Sequence[Sequence[int]],
    denominator: float,
    mode: str,
    pdf_path: Path,
) -> None:
    if mode == "density220":
        title = "Cross-treatment IBD1 sharing (>12 cM; >220 SNP/cM)"
    elif mode == "native":
        title = "Cross-treatment IBD1 sharing (>12 cM; no density filter)"
    else:
        raise ValueError(mode)

    figure, axis = plt.subplots(figsize=(8.6, 6.65))
    # All 64 cells use the same eligible-map denominator, so displaying summed
    # cM is equivalent to displaying map fraction but is easier to interpret.
    # Segment count is shown separately because a higher count can reflect
    # fragmentation rather than greater total sharing. A disclosed truncated
    # range prevents the high absolute totals from collapsing into one color.
    scale_min, scale_max = COLOR_LIMITS_CM[mode]
    totals_cm = [
        float(value) * denominator / 100.0
        for row in values
        for value in row
    ]
    if min(totals_cm) < scale_min or max(totals_cm) > scale_max:
        raise ValueError(
            f"{mode} totals fall outside declared color limits: "
            f"observed={min(totals_cm):.3f}-{max(totals_cm):.3f}, "
            f"limits={scale_min:.1f}-{scale_max:.1f}"
        )
    normalizer = Normalize(vmin=scale_min, vmax=scale_max, clip=True)
    color_map = plt.get_cmap("viridis")

    for row_index in range(8):
        for column_index in range(8):
            fraction_pct = float(values[row_index][column_index])
            total_cm = fraction_pct * denominator / 100.0
            segment_count = int(counts[row_index][column_index])
            face = color_map(normalizer(total_cm))
            axis.add_patch(
                Rectangle(
                    (column_index - 0.5, row_index - 0.5),
                    1.0,
                    1.0,
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.9,
                )
            )
            axis.text(
                column_index,
                row_index,
                f"{total_cm:,.0f}\n$n={segment_count}$",
                ha="center",
                va="center",
                fontsize=8.8,
                fontweight="semibold",
                color=annotation_color(face),
            )

    axis.set_xlim(-0.5, 7.5)
    axis.set_ylim(7.5, -0.5)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xticks(range(8))
    axis.set_yticks(range(8))
    repeated_labels = [METHOD_LABELS[method] for method in METHODS] * 2
    axis.set_xticklabels(repeated_labels, rotation=32, ha="right", fontsize=8.8)
    axis.set_yticklabels(repeated_labels, fontsize=8.8)
    axis.tick_params(which="both", length=0)

    # Strong separators show the full-UDG/non-UDG 4-by-4 blocks.
    axis.axvline(3.5, color="#222222", linewidth=2.0, zorder=4)
    axis.axhline(3.5, color="#222222", linewidth=2.0, zorder=4)
    for spine in axis.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.9)

    axis.text(
        0.25,
        1.025,
        "ORT16 full-UDG",
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.6,
        fontweight="semibold",
    )
    axis.text(
        0.75,
        1.025,
        "ORT16 non-UDG",
        transform=axis.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.6,
        fontweight="semibold",
    )
    axis.text(
        -0.155,
        0.75,
        "ORT15 full-UDG",
        transform=axis.transAxes,
        ha="center",
        va="center",
        rotation=90,
        fontsize=9.6,
        fontweight="semibold",
    )
    axis.text(
        -0.155,
        0.25,
        "ORT15 non-UDG",
        transform=axis.transAxes,
        ha="center",
        va="center",
        rotation=90,
        fontsize=9.6,
        fontweight="semibold",
    )
    scale = ScalarMappable(norm=normalizer, cmap=color_map)
    scale.set_array([])
    color_bar = figure.colorbar(
        scale,
        ax=axis,
        ticks=tuple(range(int(scale_min), int(scale_max) + 1, 100)),
        fraction=0.047,
        pad=0.045,
    )
    color_bar.set_label("IBD1 length (cM; truncated scale)", fontsize=9)
    color_bar.ax.tick_params(labelsize=8.8)

    figure.text(
        0.5,
        0.022,
        "Cell labels: summed IBD1 length (cM); n = retained segments",
        ha="center",
        va="bottom",
        fontsize=8.8,
        color="#444444",
    )

    # Count matrices are validated along with the percentages and drawn as the
    # second line in every cell.
    if len(counts) != 8 or any(len(row) != 8 for row in counts):
        raise ValueError("Count matrix is not 8-by-8")

    figure.subplots_adjust(left=0.18, right=0.91, bottom=0.13, top=0.93)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path = pdf_path.with_suffix(".png")
    metadata = {
        "Title": title,
        "Author": "ORT study",
        "Subject": "Custom summary of corrected AADR v62 ancIBD calls",
        "Creator": "plot_cross64_heatmap.py",
    }
    figure.savefig(pdf_path, metadata=metadata)
    figure.savefig(png_path, dpi=300, metadata={"Software": "Matplotlib"})
    plt.close(figure)


def write_validation(
    output_path: Path,
    summary_path: Path,
    manifest_path: Path,
    denominator: float,
    matrices: Mapping[str, Sequence[Sequence[float]]],
) -> None:
    payload = {
        "status": "PASS",
        "design": "fixed ordered 8-by-8 ORT15-by-ORT16",
        "pairs": 64,
        "resource": "AADR v62.0 exact all-SNP intersection",
        "callable_map_denominator_cM": denominator,
        "filter_semantics": {
            "length": "strict >12 cM",
            "primary_density": "strict >220 SNP/cM",
            "sensitivity_density": "none",
            "external_gap_merge": False,
        },
        "summary_sha256": sha256(summary_path),
        "resource_manifest_sha256": sha256(manifest_path),
        "ranges_pct": {
            mode: {
                "minimum": min(value for row in matrix for value in row),
                "maximum": max(value for row in matrix for value in row),
            }
            for mode, matrix in matrices.items()
        },
        "color_scales_cM": {
            mode: {"minimum": limits[0], "maximum": limits[1]}
            for mode, limits in COLOR_LIMITS_CM.items()
        },
    }
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot corrected 64-pair ancIBD treatment heatmaps."
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Corrected CROSS64_PAIR_SUMMARY.tsv",
    )
    parser.add_argument("--resource-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    denominator, matrices, counts = validate_and_build_matrices(
        args.summary, args.resource_manifest
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    primary = args.output
    draw_heatmap(
        matrices["density220"],
        counts["density220"],
        denominator,
        "density220",
        primary,
    )
    write_validation(
        args.output.with_suffix(".validation.json"),
        args.summary,
        args.resource_manifest,
        denominator,
        matrices,
    )

    print("design_validation=PASS (64 ordered pairs)")
    print(f"callable_map_denominator_cM={denominator:.6f}")
    print(f"primary={primary}")


if __name__ == "__main__":
    main()
