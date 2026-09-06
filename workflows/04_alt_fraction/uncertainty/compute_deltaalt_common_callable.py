#!/usr/bin/env python3
"""Common-callable sensitivity analysis for prespecified DeltaALT contrasts.

This is a deliberately separate sensitivity analysis.  It does not replace the
correction-specific-site primary analysis in ``compute_deltaalt_uncertainty.py``.
For each prespecified comparison, source TSVs are intersected by the exact
CHR/POS/REF/ALT key.  A site is callable only when AD_REF + AD_ALT >= 3 in every
state required by that comparison.  DeltaALT is then recomputed for every state
on that identical site set before the paired contrast is formed.

Only these comparisons are evaluated, separately for ORT15 and ORT16:

* raw non-UDG minus raw full-UDG;
* non-UDG bamRefine-5 minus non-UDG Rescale-5;
* (non-UDG Rescale-5 - raw) - (full-UDG Rescale-5 - raw); and
* (non-UDG Rescale-12 - raw) - (full-UDG Rescale-12 - raw).

Uncertainty uses the same paired leave-one-autosome-out, 5-Mb, and 10-Mb
delete-one jackknife definitions as the primary script.  Input TSVs are opened
read-only and streamed; only block-level sufficient statistics are retained.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence, Tuple

# Keep the scientific and block definitions mechanically identical to the
# primary analysis.  Deploy this file beside compute_deltaalt_uncertainty.py.
from compute_deltaalt_uncertainty import (  # type: ignore
    AUTOSOMES,
    COMMON_5MB_BLOCKS,
    COMMON_10MB_BLOCKS,
    DAMAGE,
    EXPECTED_HEADERS,
    GRCH37_AUTOSOME_LENGTHS,
    TRANSITIONS,
    ZERO_STATS,
    Stats,
    add_covered,
    add_stats,
    add_to,
    block_sort_key,
    delta,
    estimate,
    fixed_block_bounds,
    fixed_block_index,
    freeze_stats,
    jackknife_summary,
    normalize_chrom,
    output_directory,
    sha256,
    source_path,
    state_label,
    subtract_stats,
    write_csv,
)


State = Tuple[str, str, str]
Block = Tuple[str, int]
VariantKey = Tuple[int, int, str, str]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass(frozen=True)
class CommonComparison:
    comparison_id: str
    individual: str
    description: str
    states: tuple[State, ...]
    terms: tuple[tuple[State, float], ...]


def build_comparisons() -> list[CommonComparison]:
    comparisons: list[CommonComparison] = []
    for individual in ("ORT15", "ORT16"):
        def s(library: str, treatment: str) -> State:
            return individual, library, treatment

        comparisons.append(CommonComparison(
            "raw_nu_minus_raw_fu",
            individual,
            "Raw non-UDG minus raw full-UDG DeltaALT on their common callable sites",
            (s("nu", "raw"), s("fu", "raw")),
            ((s("nu", "raw"), 1.0), (s("fu", "raw"), -1.0)),
        ))
        comparisons.append(CommonComparison(
            "nu_bamrefine5_minus_rescale5",
            individual,
            "Non-UDG bamRefine-5 minus Rescale-5 DeltaALT on common callable sites",
            (s("nu", "bamrefine5"), s("nu", "rescale5")),
            ((s("nu", "bamrefine5"), 1.0), (s("nu", "rescale5"), -1.0)),
        ))
        for treatment in ("rescale5", "rescale12"):
            states = (
                s("nu", "raw"), s("nu", treatment),
                s("fu", "raw"), s("fu", treatment),
            )
            comparisons.append(CommonComparison(
                f"{treatment}_shift_nu_minus_fu",
                individual,
                f"Difference in {treatment} shift: non-UDG minus full-UDG, "
                "on four-way common callable sites",
                states,
                (
                    (s("nu", treatment), 1.0), (s("nu", "raw"), -1.0),
                    (s("fu", treatment), -1.0), (s("fu", "raw"), 1.0),
                ),
            ))
    assert len(comparisons) == 8
    assert all(set(c.states) == {state for state, _ in c.terms} for c in comparisons)
    return comparisons


COMPARISONS = build_comparisons()


@dataclass(frozen=True)
class VariantRecord:
    key: VariantKey
    chrom: str
    pos: int
    ref: str
    alt: str
    ad_ref: int
    ad_alt: int

    @property
    def depth(self) -> int:
        return self.ad_ref + self.ad_alt

    @property
    def alt_fraction(self) -> float:
        return self.ad_alt / self.depth


class CheckedReader:
    """Read and validate one sorted TSV while exposing its current record."""

    def __init__(self, state: State, path: Path):
        self.state = state
        self.path = path
        self.before = path.stat()
        self.handle = path.open("rb")
        self.digest = hashlib.sha256()
        header_raw = self.handle.readline()
        self.digest.update(header_raw)
        header = tuple(
            header_raw.decode("utf-8", errors="strict").rstrip("\r\n").split("\t")
        )
        if header not in EXPECTED_HEADERS:
            self.handle.close()
            raise ValueError(f"Unexpected header in {path}: {header!r}")
        self.header = header
        self.line_number = 1
        self.data_rows = 0
        self.unique_variant_keys = 0
        self.exact_duplicate_rows = 0
        self.seen_chroms: set[str] = set()
        self.previous_order: tuple[int, int] | None = None
        self.pending_records: list[VariantRecord] = []
        self.lookahead: tuple[VariantRecord, tuple[str, ...]] | None = None
        self.source_exhausted = False
        self.current: VariantRecord | None = None
        self.eof = False
        self.advance()

    def read_raw_record(self) -> tuple[VariantRecord, tuple[str, ...]] | None:
        """Read one raw row, enforcing only the source's CHR/POS ordering."""
        if self.source_exhausted:
            return None
        raw = self.handle.readline()
        if not raw:
            self.source_exhausted = True
            return None
        self.digest.update(raw)
        self.line_number += 1
        self.data_rows += 1
        fields = raw.decode("utf-8", errors="strict").rstrip("\r\n").split("\t")
        if len(fields) != len(self.header):
            raise ValueError(
                f"{self.path}:{self.line_number}: expected {len(self.header)} "
                f"fields, observed {len(fields)}"
            )
        try:
            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            ad_ref = int(fields[5] or 0)
            ad_alt = int(fields[6] or 0)
        except ValueError as exc:
            raise ValueError(f"{self.path}:{self.line_number}: {exc}") from exc
        if chrom not in GRCH37_AUTOSOME_LENGTHS:
            raise ValueError(f"{self.path}:{self.line_number}: non-autosomal {chrom}")
        if not 1 <= pos <= GRCH37_AUTOSOME_LENGTHS[chrom]:
            raise ValueError(f"{self.path}:{self.line_number}: position outside GRCh37")
        if ad_ref < 0 or ad_alt < 0:
            raise ValueError(f"{self.path}:{self.line_number}: negative allele count")
        order = (int(chrom), pos)
        if self.previous_order is not None and order < self.previous_order:
            raise ValueError(
                f"{self.path}:{self.line_number}: CHR/POS order decreased from "
                f"{self.previous_order[0]}:{self.previous_order[1]} to {chrom}:{pos}"
            )
        self.previous_order = order
        self.seen_chroms.add(chrom)
        ref, alt = fields[2].upper(), fields[3].upper()
        key = (int(chrom), pos, ref, alt)
        return VariantRecord(key, chrom, pos, ref, alt, ad_ref, ad_alt), tuple(fields)

    def advance(self) -> None:
        if self.eof:
            self.current = None
            return
        if self.pending_records:
            self.current = self.pending_records.pop(0)
            return

        seed = self.lookahead
        self.lookahead = None
        if seed is None:
            seed = self.read_raw_record()
        if seed is None:
            self.current = None
            self.eof = True
            return

        first_record, first_fields = seed
        coordinate = (first_record.chrom, first_record.pos)
        # REF/ALT rows at one coordinate are not guaranteed to be lexicographically
        # ordered in the inputs.  Buffer only this coordinate, validate duplicate
        # full keys, then expose a deterministic local REF/ALT ordering.
        grouped: dict[VariantKey, tuple[VariantRecord, tuple[str, ...]]] = {
            first_record.key: (first_record, first_fields)
        }
        while True:
            item = self.read_raw_record()
            if item is None:
                break
            record, row_fields = item
            if (record.chrom, record.pos) != coordinate:
                self.lookahead = item
                break
            previous = grouped.get(record.key)
            if previous is not None:
                if row_fields != previous[1]:
                    raise ValueError(
                        f"{self.path}:{self.line_number}: conflicting duplicate "
                        f"CHR/POS/REF/ALT {record.chrom}:{record.pos}:{record.ref}:{record.alt}"
                    )
                self.exact_duplicate_rows += 1
                continue
            grouped[record.key] = item

        ordered = [grouped[key][0] for key in sorted(grouped)]
        self.unique_variant_keys += len(ordered)
        self.current = ordered[0]
        self.pending_records = ordered[1:]

    def drain(self) -> None:
        while self.current is not None:
            self.advance()

    def finish(self) -> dict[str, object]:
        if not self.eof:
            self.drain()
        self.handle.close()
        if self.seen_chroms != set(AUTOSOMES):
            raise ValueError(
                f"{self.path}: chromosome set is not exactly 1..22; "
                f"missing={sorted(set(AUTOSOMES)-self.seen_chroms, key=int)}"
            )
        after = self.path.stat()
        if (self.before.st_size, self.before.st_mtime_ns) != (
            after.st_size, after.st_mtime_ns
        ):
            raise RuntimeError(f"Source changed while being read: {self.path}")
        return {
            "state": self.state,
            "source": str(self.path),
            "size_bytes": self.before.st_size,
            "mtime_ns": self.before.st_mtime_ns,
            "source_sha256": self.digest.hexdigest(),
            "data_rows": self.data_rows,
            "unique_variant_keys": self.unique_variant_keys,
            "exact_duplicate_rows_collapsed": self.exact_duplicate_rows,
        }


def new_state_accumulator() -> dict[str, object]:
    return {
        "direct": [0.0, 0, 0.0, 0, 0],
        "chrom": {},
        "5mb": {},
        "10mb": {},
    }


def record_fraction(
    accumulator: dict[str, object], chrom: str, pos: int, group: str, fraction: float
) -> None:
    direct: list[float | int] = accumulator["direct"]  # type: ignore[assignment]
    if group == "damage":
        direct[0] = float(direct[0]) + fraction
        direct[1] = int(direct[1]) + 1
    else:
        direct[2] = float(direct[2]) + fraction
        direct[3] = int(direct[3]) + 1
    direct[4] = int(direct[4]) + 1
    for scheme, block in (
        ("chrom", chrom),
        ("5mb", (chrom, fixed_block_index(chrom, pos, 5_000_000))),
        ("10mb", (chrom, fixed_block_index(chrom, pos, 10_000_000))),
    ):
        mapping: dict = accumulator[scheme]  # type: ignore[assignment]
        add_to(mapping, block, group, fraction)
        add_covered(mapping, block)


def freeze_accumulator(accumulator: Mapping[str, object]) -> dict[str, object]:
    direct = accumulator["direct"]
    return {
        "direct": (
            float(direct[0]), int(direct[1]), float(direct[2]),
            int(direct[3]), int(direct[4]),
        ),
        "chrom": freeze_stats(accumulator["chrom"]),
        "5mb": freeze_stats(accumulator["5mb"]),
        "10mb": freeze_stats(accumulator["10mb"]),
    }


def run_intersection(task: tuple[CommonComparison, str]) -> dict[str, object]:
    comparison, root_text = task
    root = Path(root_text)
    readers: list[CheckedReader] = []
    try:
        for state in comparison.states:
            readers.append(CheckedReader(state, source_path(root, state)))
    except BaseException:
        for reader in readers:
            reader.handle.close()
        raise
    accumulators = {state: new_state_accumulator() for state in comparison.states}
    diagnostics = {
        "matched_variant_keys": 0,
        "common_depth_ge3_keys": 0,
        "eligible_damage_keys": 0,
        "eligible_transversion_keys": 0,
        "excluded_other_snp_class_keys": 0,
    }
    try:
        # Sorted multi-file intersection.  A reader behind the largest current
        # key advances until it catches up; no per-site records are retained.
        while all(reader.current is not None for reader in readers):
            high = max(reader.current.key for reader in readers if reader.current is not None)
            moved = False
            for reader in readers:
                while reader.current is not None and reader.current.key < high:
                    reader.advance()
                    moved = True
                if reader.current is None:
                    break
                if reader.current.key > high:
                    high = reader.current.key
                    moved = True
            if any(reader.current is None for reader in readers):
                break
            keys = [reader.current.key for reader in readers if reader.current is not None]
            if len(set(keys)) != 1:
                # At least one reader raised the target; converge in the next pass.
                if not moved:
                    raise AssertionError("Intersection failed to make progress")
                continue

            diagnostics["matched_variant_keys"] += 1
            records = [reader.current for reader in readers]
            assert all(record is not None for record in records)
            concrete = [record for record in records if record is not None]
            if all(record.depth >= 3 for record in concrete):
                diagnostics["common_depth_ge3_keys"] += 1
                first = concrete[0]
                pair = (first.ref, first.alt)
                if pair in DAMAGE:
                    group = "damage"
                    diagnostics["eligible_damage_keys"] += 1
                elif (
                    first.ref in "ACGT" and first.alt in "ACGT"
                    and first.ref != first.alt and pair not in TRANSITIONS
                ):
                    group = "transversion"
                    diagnostics["eligible_transversion_keys"] += 1
                else:
                    group = ""
                    diagnostics["excluded_other_snp_class_keys"] += 1
                if group:
                    for state, record in zip(comparison.states, concrete):
                        fraction = record.alt_fraction
                        if not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
                            raise ValueError(f"Invalid ALT fraction for {state_label(state)}")
                        record_fraction(
                            accumulators[state], first.chrom, first.pos, group, fraction
                        )
            for reader in readers:
                reader.advance()
    finally:
        # Complete every byte stream so source hashes and chromosome gates cover
        # the entire files even if one member of the intersection ends first.
        manifests = []
        pending_error: BaseException | None = None
        for reader in readers:
            try:
                manifests.append(reader.finish())
            except BaseException as exc:  # preserve validation failures after cleanup
                if pending_error is None:
                    pending_error = exc
        if pending_error is not None:
            raise pending_error

    frozen = {state: freeze_accumulator(accumulators[state]) for state in comparison.states}
    counts = {
        (stats["direct"][1], stats["direct"][3], stats["direct"][4])
        for stats in frozen.values()
    }
    if len(counts) != 1:
        raise AssertionError("Comparison states do not have identical common-site counts")
    direct = next(iter(frozen.values()))["direct"]
    if direct[1] == 0 or direct[3] == 0:
        raise ValueError(f"{comparison.comparison_id}: a required SNP class has zero sites")
    return {
        "comparison": comparison,
        "states": frozen,
        "diagnostics": {
            **diagnostics,
            "exact_duplicate_rows_collapsed_across_input_reads": sum(
                int(row["exact_duplicate_rows_collapsed"])
                for row in manifests
            ),
        },
        "source_manifests": manifests,
    }


def block_maps(
    state_results: Mapping[State, Mapping[str, object]], scheme: str
) -> tuple[list[Block], dict[State, dict[Block, Stats]]]:
    if scheme == "loco":
        universe = [(chrom, 0) for chrom in AUTOSOMES]
        maps = {
            state: {
                (chrom, 0): result["chrom"].get(chrom, ZERO_STATS)
                for chrom in AUTOSOMES
            }
            for state, result in state_results.items()
        }
    else:
        universe = list(COMMON_5MB_BLOCKS if scheme == "5mb" else COMMON_10MB_BLOCKS)
        maps = {
            state: {block: result[scheme].get(block, ZERO_STATS) for block in universe}
            for state, result in state_results.items()
        }
    return sorted(universe, key=block_sort_key), maps


def summarize_comparison(result: Mapping[str, object]) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    comparison: CommonComparison = result["comparison"]  # type: ignore[assignment]
    state_results: Mapping[State, Mapping[str, object]] = result["states"]  # type: ignore[assignment]
    state_summary: list[dict] = []
    state_replicates: list[dict] = []
    contrast_summary: list[dict] = []
    contrast_replicates: list[dict] = []
    for scheme in ("loco", "5mb", "10mb"):
        universe, maps = block_maps(state_results, scheme)
        totals = {
            state: add_stats(maps[state][block] for block in universe)
            for state in comparison.states
        }
        for state in comparison.states:
            full = delta(totals[state])
            components = estimate(totals[state])
            assert components is not None
            deleted = []
            for block in universe:
                value = delta(subtract_stats(totals[state], maps[state][block]))
                deleted.append(value)
                state_replicates.append({
                    "scheme": scheme,
                    "individual": comparison.individual,
                    "comparison_id": comparison.comparison_id,
                    "library": state[1],
                    "treatment": state[2],
                    "omitted_chrom": block[0],
                    "omitted_block_index_0based": block[1],
                    "delete_one_deltaalt_pp": value,
                })
            state_summary.append({
                "scheme": scheme,
                "individual": comparison.individual,
                "comparison_id": comparison.comparison_id,
                "library": state[1],
                "treatment": state[2],
                "damage_alt_fraction_sum": totals[state][0],
                "damage_n_common_dp3": totals[state][1],
                "transversion_alt_fraction_sum": totals[state][2],
                "transversion_n_common_dp3": totals[state][3],
                "common_target_sites": totals[state][4],
                "damage_alt_pct": components[0],
                "transversion_alt_pct": components[1],
                "estimate_unit": "percentage_points",
                **jackknife_summary(full, deleted),
            })

        full_contrast = math.fsum(
            coefficient * delta(totals[state])
            for state, coefficient in comparison.terms
        )
        deleted_contrasts = []
        for block in universe:
            value = math.fsum(
                coefficient * delta(subtract_stats(totals[state], maps[state][block]))
                for state, coefficient in comparison.terms
            )
            deleted_contrasts.append(value)
            contrast_replicates.append({
                "scheme": scheme,
                "individual": comparison.individual,
                "comparison_id": comparison.comparison_id,
                "omitted_chrom": block[0],
                "omitted_block_index_0based": block[1],
                "delete_one_contrast_pp": value,
            })
        contrast_summary.append({
            "scheme": scheme,
            "individual": comparison.individual,
            "comparison_id": comparison.comparison_id,
            "description": comparison.description,
            "common_callable_definition": (
                "exact CHR/POS/REF/ALT intersection and AD_REF+AD_ALT>=3 in every term"
            ),
            "terms_json": json.dumps([
                {"state": state_label(state), "coefficient": coefficient}
                for state, coefficient in comparison.terms
            ], separators=(",", ":")),
            "estimate_unit": "percentage_points",
            **jackknife_summary(full_contrast, deleted_contrasts),
        })
    return state_summary, state_replicates, contrast_summary, contrast_replicates


def sufficient_rows(result: Mapping[str, object], scheme: str) -> list[dict]:
    comparison: CommonComparison = result["comparison"]  # type: ignore[assignment]
    state_results: Mapping[State, Mapping[str, object]] = result["states"]  # type: ignore[assignment]
    rows: list[dict] = []
    universe, maps = block_maps(state_results, "loco" if scheme == "chrom" else scheme)
    size = None if scheme == "chrom" else (5_000_000 if scheme == "5mb" else 10_000_000)
    for state in comparison.states:
        for chrom, index in universe:
            stats = maps[state][(chrom, index)]
            if size is None:
                block_id, start, end = f"chr{chrom}", "", ""
            else:
                start, end = fixed_block_bounds(chrom, index, size)
                block_id = f"chr{chrom}:{start}-{end}"
            components = estimate(stats)
            rows.append({
                "individual": comparison.individual,
                "comparison_id": comparison.comparison_id,
                "library": state[1],
                "treatment": state[2],
                "scheme": scheme,
                "chrom": chrom,
                "block_id": block_id,
                "block_index_0based": "" if size is None else index,
                "block_start_1based": start,
                "block_end_inclusive": end,
                "damage_alt_fraction_sum": stats[0],
                "damage_n_common_dp3": stats[1],
                "transversion_alt_fraction_sum": stats[2],
                "transversion_n_common_dp3": stats[3],
                "common_target_sites": stats[4],
                "deltaalt_pp": "" if components is None else components[2],
            })
    return rows


def check_internal(result: Mapping[str, object]) -> None:
    comparison: CommonComparison = result["comparison"]  # type: ignore[assignment]
    state_results: Mapping[State, Mapping[str, object]] = result["states"]  # type: ignore[assignment]
    direct_counts = set()
    for state, state_result in state_results.items():
        direct: Stats = state_result["direct"]  # type: ignore[assignment]
        direct_counts.add((direct[1], direct[3], direct[4]))
        for scheme, universe in (
            ("chrom", [(chrom, 0) for chrom in AUTOSOMES]),
            ("5mb", COMMON_5MB_BLOCKS),
            ("10mb", COMMON_10MB_BLOCKS),
        ):
            if scheme == "chrom":
                observed = state_result[scheme]
                combined = add_stats(observed.get(block[0], ZERO_STATS) for block in universe)
            else:
                observed = state_result[scheme]
                combined = add_stats(observed.get(block, ZERO_STATS) for block in universe)
            if (combined[1], combined[3], combined[4]) != (direct[1], direct[3], direct[4]):
                raise AssertionError(
                    f"{comparison.comparison_id} {state_label(state)} {scheme} count mismatch"
                )
            tolerance = max(1e-8, 2e-12 * max(abs(direct[0]), abs(direct[2]), 1.0))
            if not math.isclose(combined[0], direct[0], rel_tol=0, abs_tol=tolerance):
                raise AssertionError("Damage sums do not reconstruct")
            if not math.isclose(combined[2], direct[2], rel_tol=0, abs_tol=tolerance):
                raise AssertionError("Transversion sums do not reconstruct")
        # Because the site set is common by construction, every physical block
        # must contain the same class counts for every state in the comparison.
        for scheme in ("chrom", "5mb", "10mb"):
            universe, maps = block_maps(
                state_results, "loco" if scheme == "chrom" else scheme
            )
            for block in universe:
                counts = {
                    (maps[s][block][1], maps[s][block][3], maps[s][block][4])
                    for s in comparison.states
                }
                if len(counts) != 1:
                    raise AssertionError(
                        f"{comparison.comparison_id} {scheme} {block}: "
                        "common-site block counts differ by state"
                    )
    if len(direct_counts) != 1:
        raise AssertionError(f"{comparison.comparison_id}: common-site counts differ by state")


def block_class_audit_rows(result: Mapping[str, object]) -> list[dict]:
    """Report class support once per comparison and physical block."""
    comparison: CommonComparison = result["comparison"]  # type: ignore[assignment]
    state_results: Mapping[State, Mapping[str, object]] = result["states"]  # type: ignore[assignment]
    rows: list[dict] = []
    reference_state = comparison.states[0]
    for scheme in ("loco", "5mb", "10mb"):
        universe, maps = block_maps(state_results, scheme)
        size = None if scheme == "loco" else (5_000_000 if scheme == "5mb" else 10_000_000)
        total = add_stats(maps[reference_state][block] for block in universe)
        for chrom, index in universe:
            block = (chrom, index)
            stats = maps[reference_state][block]
            if size is None:
                block_id, start, end = f"chr{chrom}", 1, GRCH37_AUTOSOME_LENGTHS[chrom]
            else:
                start, end = fixed_block_bounds(chrom, index, size)
                block_id = f"chr{chrom}:{start}-{end}"
            rows.append({
                "individual": comparison.individual,
                "comparison_id": comparison.comparison_id,
                "scheme": scheme,
                "chrom": chrom,
                "block_id": block_id,
                "block_index_0based": "" if size is None else index,
                "block_start_1based": start,
                "block_end_inclusive": end,
                "damage_n_common_dp3": stats[1],
                "transversion_n_common_dp3": stats[3],
                "common_target_sites": stats[4],
                "block_has_damage_sites": stats[1] > 0,
                "block_has_transversion_sites": stats[3] > 0,
                "block_has_both_classes": stats[1] > 0 and stats[3] > 0,
                "delete_one_has_both_classes": (
                    total[1] - stats[1] > 0 and total[3] - stats[3] > 0
                ),
            })
    return rows


def run_self_tests() -> None:
    assert len(COMPARISONS) == 8
    assert len({(c.individual, c.comparison_id) for c in COMPARISONS}) == 8
    assert {len(c.states) for c in COMPARISONS} == {2, 4}
    probe: Stats = (3.0, 2, 1.0, 2, 4)
    assert delta(probe) == 100.0
    for comparison in COMPARISONS:
        assert math.isclose(math.fsum(coef for _, coef in comparison.terms), 0.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="New or empty output directory")
    parser.add_argument(
        "--source-root", type=Path, required=True,
        help="External GLIMPSE analysis root containing the required per-site TSVs",
    )
    parser.add_argument("--workers", type=int, default=8, help="Parallel intersections (1-8)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_self_tests()
    if not 1 <= args.workers <= 8:
        raise ValueError("--workers must be between 1 and 8")
    unique_states = {state for comparison in COMPARISONS for state in comparison.states}
    missing = [source_path(args.source_root, state) for state in unique_states
               if not source_path(args.source_root, state).is_file()]
    if missing:
        raise FileNotFoundError("Missing source TSV(s):\n" + "\n".join(map(str, missing)))
    output_directory(args.output_dir)
    started = utc_now()
    results: dict[tuple[str, str], Mapping[str, object]] = {}
    tasks = [(comparison, str(args.source_root)) for comparison in COMPARISONS]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_intersection, task): task[0] for task in tasks}
        for future in as_completed(futures):
            comparison = futures[future]
            result = future.result()
            key = (comparison.individual, comparison.comparison_id)
            results[key] = result
            print(f"completed {comparison.individual}:{comparison.comparison_id}",
                  file=sys.stderr, flush=True)
    if len(results) != 8:
        raise AssertionError("Did not receive all eight comparison results")

    ordered = [results[(c.individual, c.comparison_id)] for c in COMPARISONS]
    for result in ordered:
        check_internal(result)

    block_audit: list[dict] = []
    for result in ordered:
        block_audit.extend(block_class_audit_rows(result))
    if not all(bool(row["delete_one_has_both_classes"]) for row in block_audit):
        raise ValueError(
            "At least one delete-one replicate has an undefined SNP-class denominator"
        )
    write_csv(args.output_dir / "common_callable_block_class_audit.csv", block_audit)

    state_summary: list[dict] = []
    state_rep: list[dict] = []
    contrast_summary: list[dict] = []
    contrast_rep: list[dict] = []
    for result in ordered:
        a, b, c, d = summarize_comparison(result)
        state_summary.extend(a)
        state_rep.extend(b)
        contrast_summary.extend(c)
        contrast_rep.extend(d)
    write_csv(args.output_dir / "common_callable_state_summary.csv", state_summary)
    write_csv(args.output_dir / "common_callable_state_delete_one.csv", state_rep)
    write_csv(args.output_dir / "common_callable_contrast_summary.csv", contrast_summary)
    write_csv(args.output_dir / "common_callable_contrast_delete_one.csv", contrast_rep)

    for scheme in ("chrom", "5mb", "10mb"):
        rows = []
        for result in ordered:
            rows.extend(sufficient_rows(result, scheme))
        write_csv(args.output_dir / f"common_callable_sufficient_stats_{scheme}.csv", rows)

    diagnostics = []
    source_observations: dict[str, list[Mapping[str, object]]] = {}
    for result in ordered:
        comparison: CommonComparison = result["comparison"]  # type: ignore[assignment]
        diagnostics.append({
            "individual": comparison.individual,
            "comparison_id": comparison.comparison_id,
            "n_states_intersected": len(comparison.states),
            **result["diagnostics"],
        })
        for manifest in result["source_manifests"]:
            source_observations.setdefault(str(manifest["source"]), []).append(manifest)
    write_csv(args.output_dir / "common_callable_diagnostics.csv", diagnostics)

    source_rows = []
    for path, observations in sorted(source_observations.items()):
        states = {tuple(row["state"]) for row in observations}
        digests = {str(row["source_sha256"]) for row in observations}
        sizes = {int(row["size_bytes"]) for row in observations}
        mtimes = {int(row["mtime_ns"]) for row in observations}
        data_rows = {int(row["data_rows"]) for row in observations}
        unique_keys = {int(row["unique_variant_keys"]) for row in observations}
        duplicates = {
            int(row["exact_duplicate_rows_collapsed"])
            for row in observations
        }
        if any(len(values) != 1 for values in (
            states, digests, sizes, mtimes, data_rows, unique_keys, duplicates
        )):
            raise RuntimeError(f"Repeated source observations disagree: {path}")
        state = next(iter(states))
        source_rows.append({
            "individual": state[0],
            "library": state[1],
            "treatment": state[2],
            "state": state_label(state),
            "source_file": Path(path).name,
            "size_bytes": next(iter(sizes)),
            "source_sha256": next(iter(digests)),
            "data_rows": next(iter(data_rows)),
            "unique_variant_keys": next(iter(unique_keys)),
            "exact_duplicate_rows_collapsed": next(iter(duplicates)),
            "times_streamed": len(observations),
            "source_modified": "no",
        })
    expected_states = {state for comparison in COMPARISONS for state in comparison.states}
    manifested_states = {
        (str(row["individual"]), str(row["library"]), str(row["treatment"]))
        for row in source_rows
    }
    if len(source_rows) != 14 or manifested_states != expected_states:
        raise RuntimeError("Common-callable source manifest is not the expected 14 states")
    write_csv(args.output_dir / "common_callable_source_manifest.csv", source_rows)

    duplicate_total = sum(
        int(row["exact_duplicate_rows_collapsed"])
        for row in source_rows
    )
    blocks_without_both = sum(
        not bool(row["block_has_both_classes"]) for row in block_audit
    )

    validation_lines = [
        "# Common-callable DeltaALT sensitivity validation\n",
        "- Scope: eight prespecified contrasts only (four per individual).",
        "- Exact variant-key intersection: CHR/POS/REF/ALT.",
        "- Common callability: AD_REF + AD_ALT >= 3 in every comparison state.",
        "- Common target-site counts: identical across all states within each comparison.",
        f"- Exact-identical duplicate rows collapsed within coordinate groups: {duplicate_total}.",
        "- Conflicting duplicate CHR/POS/REF/ALT keys: rejected.",
        "- Uncertainty: paired LOCO, 5-Mb, and 10-Mb delete-one jackknives.",
        "- Physical block grids and t-based intervals: identical to the primary script.",
        f"- Comparison-blocks without both SNP classes: {blocks_without_both}; recorded in",
        "  common_callable_block_class_audit.csv. Every delete-one denominator was defined.",
        "- Entire-source validation: exact schema, monotonic CHR/POS order with locally",
        "  sorted same-position REF/ALT keys, chromosome set 1-22,",
        "  nonnegative allele counts, full-stream SHA-256, and unchanged size/mtime.",
        "- Source files modified: **no**.",
        "- This sensitivity analysis does not replace the correction-specific-site primary analysis.\n",
    ]
    validation_path = args.output_dir / "VALIDATION.md"
    tmp_validation = validation_path.with_suffix(".md.tmp")
    with tmp_validation.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(validation_lines))
    os.replace(tmp_validation, validation_path)

    script = Path(__file__).resolve()
    manifest = {
        "status": "complete",
        "started_utc": started,
        "finished_utc": utc_now(),
        "platform": platform.platform(),
        "python": sys.version,
        "script": script.name,
        "script_sha256": sha256(script),
        "primary_definition_script": Path(
            sys.modules["compute_deltaalt_uncertainty"].__file__
        ).name,
        "primary_definition_script_sha256": sha256(
            Path(sys.modules["compute_deltaalt_uncertainty"].__file__).resolve()
        ),
        "source_root": "external; paths intentionally not serialized",
        "workers": args.workers,
        "comparisons": 8,
        "unique_source_files": len(source_rows),
        "source_sha256_by_state": {
            str(row["state"]): str(row["source_sha256"]) for row in source_rows
        },
        "common_callable_definition": (
            "exact CHR/POS/REF/ALT intersection; AD_REF+AD_ALT>=3 in every state "
            "required by the comparison"
        ),
        "block_schemes": ["leave-one-autosome-out", "5mb", "10mb"],
        "paired_deletion": "same genomic block omitted from every contrast term",
        "duplicate_policy": (
            "group rows by monotonic CHR/POS, locally sort REF/ALT keys, collapse "
            "exact-identical duplicate rows, and reject conflicting duplicate keys"
        ),
        "exact_duplicate_rows_collapsed": duplicate_total,
        "block_class_audit": {
            "comparison_blocks": len(block_audit),
            "blocks_without_both_classes": blocks_without_both,
            "all_delete_one_denominators_defined": True,
        },
        "source_files_modified": False,
        "internal_validation_pass": True,
    }
    tmp_manifest = args.output_dir / "run_manifest.json.tmp"
    with tmp_manifest.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_manifest, args.output_dir / "run_manifest.json")
    print(f"complete: {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
