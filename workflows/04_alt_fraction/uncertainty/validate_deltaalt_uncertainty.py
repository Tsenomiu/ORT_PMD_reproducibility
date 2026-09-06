#!/usr/bin/env python3
"""Validate compact Round45 DeltaALT uncertainty release tables."""

from __future__ import annotations

import csv
import hashlib
import math
import re
from collections import defaultdict
from pathlib import Path


RELEASE_ROOT = Path(__file__).resolve().parents[3]
SUMMARY = RELEASE_ROOT / "data/summary/damage/deltaalt_uncertainty"
CANONICAL = RELEASE_ROOT / "data/summary/damage/alt_fraction_recomputed_dp3_adsum.csv"
SCHEMES = {"loco": 22, "5mb": 577, "10mb": 289}
CANONICAL_SHA256 = "8e4d4699b18435d78c85c79a2fa2dece69dbb9a5580f24953b0e6dd04b7f1b73"
PRIMARY_CONTRAST_IDS = {
    "raw_nu_minus_raw_fu",
    "fu_trim5_minus_raw", "fu_rescale5_minus_raw",
    "fu_rescale12_minus_raw", "fu_bamrefine5_minus_raw",
    "nu_trim5_minus_raw", "nu_trim10_minus_raw",
    "nu_rescale5_minus_raw", "nu_rescale10_minus_raw",
    "nu_rescale12_minus_raw", "nu_bamrefine5_minus_raw",
    "nu_bamrefine10_minus_raw", "nu_bamrefine5_minus_rescale5",
    "rescale5_shift_nu_minus_fu", "rescale12_shift_nu_minus_fu",
}
COMMON_COMPARISON_IDS = {
    "raw_nu_minus_raw_fu", "nu_bamrefine5_minus_rescale5",
    "rescale5_shift_nu_minus_fu", "rescale12_shift_nu_minus_fu",
}
ZERO_CROSSING_LOCO_STATES = {
    ("ORT15", "fu", "bamrefine5"),
    ("ORT16", "nu", "trim10"),
    ("ORT16", "nu", "bamrefine10"),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise AssertionError(f"No rows in {path}")
    return rows


def row_index(rows: list[dict[str, str]], *keys: str) -> dict[tuple[str, ...], dict[str, str]]:
    index: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[column] for column in keys)
        if key in index:
            raise AssertionError(f"Duplicate key {key}")
        index[key] = row
    return index


def close(actual: str, expected: float, tolerance: float = 5e-12) -> None:
    if not math.isclose(float(actual), expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"Expected {expected}, observed {actual}")


def interval_crosses_zero(row: dict[str, str]) -> bool:
    return float(row["ci95_lo_pp"]) <= 0.0 <= float(row["ci95_hi_pp"])


def main() -> None:
    if hashlib.sha256(CANONICAL.read_bytes()).hexdigest() != CANONICAL_SHA256:
        raise AssertionError("Canonical 26-state DeltaALT CSV SHA-256 changed")
    canonical = row_index(read_csv(CANONICAL), "individual", "library", "treatment")
    if len(canonical) != 26:
        raise AssertionError(f"Expected 26 canonical states, observed {len(canonical)}")

    validation = row_index(
        read_csv(SUMMARY / "canonical_validation.csv"),
        "individual", "library", "treatment",
    )
    if set(validation) != set(canonical):
        raise AssertionError("Canonical validation state set differs from the 26-row source table")
    for key, row in validation.items():
        if row["validation_pass"] != "yes" or row["source_basename_match"] != "yes":
            raise AssertionError(f"Canonical validation failed for {key}")
        if int(row["malformed_rows"]) != 0:
            raise AssertionError(f"Malformed source rows reported for {key}")
        close(row["computed_deltaalt_pp"], float(canonical[key]["ts_minus_tv_pp"]), 5.1e-7)

    with (SUMMARY / "source_checksums.tsv").open(newline="", encoding="utf-8") as handle:
        checksum_rows = list(csv.DictReader(handle, delimiter="\t"))
    checksum_index = row_index(checksum_rows, "individual", "library", "treatment")
    if set(checksum_index) != set(canonical):
        raise AssertionError("Source-checksum state set differs from the 26 canonical states")
    for key, row in checksum_index.items():
        if not re.fullmatch(r"[0-9a-f]{64}", row["source_sha256"]):
            raise AssertionError(f"Invalid source SHA-256 for {key}")
        if min(int(row[field]) for field in ("size_bytes", "data_rows", "covered_sites")) <= 0:
            raise AssertionError(f"Invalid source counts for {key}")

    state_ci = row_index(
        read_csv(SUMMARY / "state_ci_loco.csv"),
        "individual", "library", "treatment",
    )
    if set(state_ci) != set(canonical):
        raise AssertionError("LOCO interval state set differs from the canonical state set")
    crossing = set()
    for key, row in state_ci.items():
        close(row["deltaalt_point_pp"], float(canonical[key]["ts_minus_tv_pp"]), 5.1e-7)
        crosses = interval_crosses_zero(row)
        if crosses:
            crossing.add(key)
        if (row["interval_excludes_zero"] == "yes") == crosses:
            raise AssertionError(f"Interval zero flag inconsistent for {key}")
    if crossing != ZERO_CROSSING_LOCO_STATES:
        raise AssertionError(f"Unexpected LOCO zero-crossing states: {sorted(crossing)}")

    state_summary = row_index(
        read_csv(SUMMARY / "primary_state_summary.csv"),
        "scheme", "individual", "library", "treatment",
    )
    if len(state_summary) != 78:
        raise AssertionError(f"Expected 78 primary state rows, observed {len(state_summary)}")
    per_state: dict[tuple[str, str, str], set[float]] = defaultdict(set)
    for (scheme, individual, library, treatment), row in state_summary.items():
        if scheme not in SCHEMES or int(row["n_blocks"]) != SCHEMES[scheme]:
            raise AssertionError(f"Unexpected block scheme/count for {(scheme, individual, library, treatment)}")
        per_state[(individual, library, treatment)].add(float(row["full_estimate_pp"]))
    if set(per_state) != set(canonical) or any(len(values) != 1 for values in per_state.values()):
        raise AssertionError("Primary point estimates are not invariant across block schemes")
    for key, values in per_state.items():
        close(str(next(iter(values))), float(canonical[key]["ts_minus_tv_pp"]), 5.1e-7)

    primary = row_index(
        read_csv(SUMMARY / "primary_paired_contrasts.csv"),
        "scheme", "individual", "contrast_id",
    )
    if len(primary) != 90:
        raise AssertionError(f"Expected 90 primary contrast rows, observed {len(primary)}")
    if {key[2] for key in primary} != PRIMARY_CONTRAST_IDS:
        raise AssertionError("Unexpected primary contrast identifier set")
    primary_invariance: dict[tuple[str, str], set[tuple[int, str]]] = defaultdict(set)
    for (scheme, individual, contrast_id), row in primary.items():
        if scheme not in SCHEMES or int(row["n_blocks"]) != SCHEMES[scheme]:
            raise AssertionError(f"Unexpected block scheme/count for {(scheme, individual, contrast_id)}")
        if contrast_id.endswith("_minus_raw"):
            if float(row["full_estimate_pp"]) >= 0 or float(row["ci95_hi_pp"]) >= 0:
                raise AssertionError(f"Correction-versus-raw interval not wholly below zero: {(scheme, individual, contrast_id)}")
        sign = 1 if float(row["full_estimate_pp"]) > 0 else -1
        primary_invariance[(individual, contrast_id)].add((sign, row["ci_excludes_zero"]))
    if any(len(values) != 1 for values in primary_invariance.values()):
        raise AssertionError("Primary contrast sign/zero-exclusion status changes by block scheme")
    for scheme in SCHEMES:
        correction_rows = [
            row for (row_scheme, _individual, contrast_id), row in primary.items()
            if row_scheme == scheme and contrast_id.endswith("_minus_raw")
        ]
        if len(correction_rows) != 22:
            raise AssertionError(f"Expected 22 correction-minus-raw contrasts for {scheme}")
    raw15 = primary[("loco", "ORT15", "raw_nu_minus_raw_fu")]
    raw16 = primary[("loco", "ORT16", "raw_nu_minus_raw_fu")]
    for row, expected in ((raw15, (0.381700786159, 0.287402361712, 0.480666006617)),
                          (raw16, (0.300258220287, 0.19203243564, 0.40918643743))):
        close(row["full_estimate_pp"], expected[0])
        close(row["ci95_lo_pp"], expected[1])
        close(row["ci95_hi_pp"], expected[2])

    common = row_index(
        read_csv(SUMMARY / "common_callable_sensitivity.csv"),
        "scheme", "individual", "comparison_id",
    )
    if len(common) != 24:
        raise AssertionError(f"Expected 24 common-callable contrast rows, observed {len(common)}")
    if {key[2] for key in common} != COMMON_COMPARISON_IDS:
        raise AssertionError("Unexpected common-callable comparison identifier set")
    common_invariance: dict[tuple[str, str], set[tuple[int, str]]] = defaultdict(set)
    for (scheme, individual, comparison_id), row in common.items():
        if scheme not in SCHEMES or int(row["n_blocks"]) != SCHEMES[scheme]:
            raise AssertionError(f"Unexpected common-callable block count: {(scheme, individual, comparison_id)}")
        sign = 1 if float(row["full_estimate_pp"]) > 0 else -1
        common_invariance[(individual, comparison_id)].add((sign, row["ci_excludes_zero"]))
    if any(len(values) != 1 for values in common_invariance.values()):
        raise AssertionError("Common-callable sign/zero-exclusion status changes by block scheme")
    for scheme, n_blocks in SCHEMES.items():
        ort15_raw = common[(scheme, "ORT15", "raw_nu_minus_raw_fu")]
        ort16_raw = common[(scheme, "ORT16", "raw_nu_minus_raw_fu")]
        if not interval_crosses_zero(ort15_raw) or interval_crosses_zero(ort16_raw):
            raise AssertionError(f"Unexpected common-callable raw contrast at {scheme}")
        primary_r12 = primary[(scheme, "ORT16", "rescale12_shift_nu_minus_fu")]
        common_r12 = common[(scheme, "ORT16", "rescale12_shift_nu_minus_fu")]
        if float(primary_r12["full_estimate_pp"]) >= 0 or float(common_r12["full_estimate_pp"]) <= 0:
            raise AssertionError(f"ORT16 Rescale-12 sensitivity did not reverse direction at {scheme}")

    common_states = read_csv(SUMMARY / "common_callable_state_summary.csv")
    diagnostics = row_index(
        read_csv(SUMMARY / "common_callable_diagnostics.csv"),
        "individual", "comparison_id",
    )
    if len(common_states) != 72 or len(diagnostics) != 8:
        raise AssertionError("Unexpected common-callable state/diagnostic row count")
    if any(int(row["common_depth_ge3_keys"]) <= 0 for row in diagnostics.values()):
        raise AssertionError("A common-callable comparison contains no eligible sites")
    expected_diagnostics = {
        ("ORT15", "raw_nu_minus_raw_fu"): (2, 19429, 15033),
        ("ORT15", "nu_bamrefine5_minus_rescale5"): (2, 164795, 132373),
        ("ORT15", "rescale5_shift_nu_minus_fu"): (4, 18639, 14892),
        ("ORT15", "rescale12_shift_nu_minus_fu"): (4, 18161, 14805),
        ("ORT16", "raw_nu_minus_raw_fu"): (2, 79042, 62293),
        ("ORT16", "nu_bamrefine5_minus_rescale5"): (2, 291091, 237808),
        ("ORT16", "rescale5_shift_nu_minus_fu"): (4, 76166, 61761),
        ("ORT16", "rescale12_shift_nu_minus_fu"): (4, 74556, 61541),
    }
    for key, expected in expected_diagnostics.items():
        row = diagnostics[key]
        observed = tuple(int(row[field]) for field in (
            "n_states_intersected", "eligible_damage_keys", "eligible_transversion_keys"
        ))
        if observed != expected:
            raise AssertionError(f"Unexpected common-callable diagnostics for {key}: {observed}")

    print("PASS: DeltaALT uncertainty compact release (26 states; 90 primary and 24 sensitivity contrasts)")


if __name__ == "__main__":
    main()
