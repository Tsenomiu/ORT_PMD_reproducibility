#!/usr/bin/env python3
"""Recompute DeltaALT and estimate genomic-block jackknife uncertainty.

This program is intentionally read-only with respect to the 26 primary TEST-server
TSVs.  It streams each TSV once, retains only block-level sufficient statistics,
validates the full-scope results against the canonical 26-row CSV, and then writes
auditable state and paired-contrast jackknife outputs.

Scientific definition
---------------------
At sites with AD_REF + AD_ALT >= 3, the site ALT fraction is
AD_ALT / (AD_REF + AD_ALT).  DeltaALT (percentage points) is 100 times the
difference between the mean fraction at C>T/G>A SNPs and the mean at
transversion SNPs.  Each correction state uses its own eligible-site set.

Primary uncertainty design
--------------------------
* leave-one-autosome-out (LOCO; 22 paired chromosome blocks), and
* paired delete-one 5-Mb and 10-Mb autosomal genomic-block jackknives.

For a paired contrast, the same genomic block is omitted from every term.  All
states and contrasts use one full GRCh37 physical-block universe per scheme;
sites are not intersected between correction states.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence, Tuple


RELEASE_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CANONICAL = (
    RELEASE_ROOT / "data/summary/damage/alt_fraction_recomputed_dp3_adsum.csv"
)
AUTOSOMES = tuple(str(x) for x in range(1, 23))
# GRCh37/hg19 autosome lengths.  The input SNP coordinates and canonical 1240k
# resources are on this build.  These lengths define a common block grid across
# every correction state.
GRCH37_AUTOSOME_LENGTHS = {
    "1": 249_250_621, "2": 243_199_373, "3": 198_022_430,
    "4": 191_154_276, "5": 180_915_260, "6": 171_115_067,
    "7": 159_138_663, "8": 146_364_022, "9": 141_213_431,
    "10": 135_534_747, "11": 135_006_516, "12": 133_851_895,
    "13": 115_169_878, "14": 107_349_540, "15": 102_531_392,
    "16": 90_354_753, "17": 81_195_210, "18": 78_077_248,
    "19": 59_128_983, "20": 63_025_520, "21": 48_129_895,
    "22": 51_304_566,
}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
DAMAGE = {("C", "T"), ("G", "A")}
ROUND45_HEADER = (
    "CHR", "POS", "REF", "ALT", "DP", "REF_COUNT", "ALT_COUNT",
    "ALT_FRAC", "CLASS",
)
PUBLIC_ALT_COUNT_HEADER = (
    "CHROM", "POS", "REF", "ALT", "DP", "AD_REF", "AD_ALT",
)
EXPECTED_HEADERS = (ROUND45_HEADER, PUBLIC_ALT_COUNT_HEADER)
State = Tuple[str, str, str]
Block = Tuple[str, int]
# Sufficient statistics: damage sum/n, transversion sum/n, and covered-site n.
Stats = Tuple[float, int, float, int, int]
ZERO_STATS: Stats = (0.0, 0, 0.0, 0, 0)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def state_label(state: State) -> str:
    return ":".join(state)


def source_path(root: Path, state: State) -> Path:
    """Return the exact source mapping used by tmp/recompute_alt_fraction.py."""
    individual, library, treatment = state
    base = root / "alt_mpileup/alt_frac_from_bcf/tables"
    rescale = root / "rescaled/alt_mpileup/alt_frac_from_bcf/tables"
    rescale_fu = root / "rescaled_fu/alt_mpileup/alt_frac_from_bcf/tables"
    stem = f"{individual}_{library}.1"
    if treatment == "raw":
        return base / f"{stem}.altfrac.q30Q20.tsv"
    if treatment == "trim5":
        return base / f"{stem}.trim5.altfrac.q30Q20.tsv"
    if treatment == "trim10":
        return base / f"{stem}.trim10.altfrac.q30Q20.tsv"
    if treatment == "bamrefine5":
        return base / f"{stem}.bamrefine5.altfrac.q30Q20.tsv"
    if treatment == "bamrefine10":
        return base / f"{stem}.bamrefine10.altfrac.q30Q20.tsv"
    if treatment == "rescale5":
        return rescale / f"{stem}.rescaled5.altfrac.q30Q20.tsv"
    if treatment == "rescale10":
        return rescale / f"{stem}.rescaled10.altfrac.q30Q20.tsv"
    if treatment == "rescale12" and library == "fu":
        return rescale_fu / f"{stem}.rescaled.altfrac.q30Q20.tsv"
    if treatment == "rescale12" and library == "nu":
        return base / f"{stem}.rescaled.altfrac.q30Q20.tsv"
    raise ValueError(f"No source mapping for {state!r}")


def build_states() -> list[State]:
    rows: list[State] = []
    for individual in ("ORT15", "ORT16"):
        for treatment in ("raw", "trim5", "rescale5", "rescale12", "bamrefine5"):
            rows.append((individual, "fu", treatment))
    for individual in ("ORT15", "ORT16"):
        for treatment in (
            "raw", "trim5", "trim10", "rescale5", "rescale10",
            "rescale12", "bamrefine5", "bamrefine10",
        ):
            rows.append((individual, "nu", treatment))
    assert len(rows) == 26 and len(set(rows)) == 26
    return rows


STATES = build_states()


def normalize_chrom(text: str) -> str:
    chrom = text.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    upper = chrom.upper()
    aliases = {"23": "X", "24": "Y", "25": "MT", "M": "MT"}
    upper = aliases.get(upper, upper)
    if upper in {"X", "Y", "MT"}:
        return upper
    try:
        value = int(upper)
    except ValueError as exc:
        raise ValueError(f"Unsupported chromosome label {text!r}") from exc
    if not 1 <= value <= 22:
        raise ValueError(f"Unsupported chromosome label {text!r}")
    return str(value)


def fixed_block_index(chrom: str, pos: int, size: int) -> int:
    """Map a GRCh37 position to a nominal fixed block, merging a short tail.

    A terminal remainder shorter than half the target width is merged into the
    preceding block.  Thus no terminal block is less than 2.5 Mb in the 5-Mb
    scheme or less than 5 Mb in the 10-Mb scheme.  The rule and reference lengths
    are fixed globally, not inferred separately from correction-specific sites.
    """
    length = GRCH37_AUTOSOME_LENGTHS[chrom]
    if not 1 <= pos <= length:
        raise ValueError(f"Position {chrom}:{pos} is outside GRCh37")
    index = (pos - 1) // size
    full_blocks, remainder = divmod(length, size)
    if remainder and remainder < size / 2 and index == full_blocks:
        index -= 1
    return index


def fixed_block_bounds(chrom: str, index: int, size: int) -> tuple[int, int]:
    length = GRCH37_AUTOSOME_LENGTHS[chrom]
    full_blocks, remainder = divmod(length, size)
    start = index * size + 1
    if remainder and remainder < size / 2 and index == full_blocks - 1:
        end = length
    else:
        end = min((index + 1) * size, length)
    return start, end


def common_fixed_blocks(size: int) -> tuple[Block, ...]:
    blocks: list[Block] = []
    for chrom in AUTOSOMES:
        length = GRCH37_AUTOSOME_LENGTHS[chrom]
        indices = {fixed_block_index(chrom, pos, size) for pos in range(1, length + 1, size)}
        indices.add(fixed_block_index(chrom, length, size))
        blocks.extend((chrom, index) for index in sorted(indices))
    return tuple(blocks)


COMMON_5MB_BLOCKS = common_fixed_blocks(5_000_000)
COMMON_10MB_BLOCKS = common_fixed_blocks(10_000_000)


def add_to(mapping: dict, key: object, group: str, fraction: float) -> None:
    values = mapping.get(key)
    if values is None:
        values = [0.0, 0, 0.0, 0, 0]
        mapping[key] = values
    if group == "damage":
        values[0] += fraction
        values[1] += 1
    else:
        values[2] += fraction
        values[3] += 1


def add_covered(mapping: dict, key: object) -> None:
    values = mapping.get(key)
    if values is None:
        values = [0.0, 0, 0.0, 0, 0]
        mapping[key] = values
    values[4] += 1


def freeze_stats(mapping: Mapping[object, Sequence[float | int]]) -> dict:
    return {
        key: (
            float(value[0]), int(value[1]), float(value[2]), int(value[3]), int(value[4])
        )
        for key, value in mapping.items()
    }


def stream_one(task: tuple[State, str]) -> dict[str, object]:
    """Worker: stream one TSV and return only compact sufficient statistics."""
    state, path_text = task
    path = Path(path_text)
    before = path.stat()
    chrom_stats: dict[str, list[float | int]] = {}
    block5_stats: dict[Block, list[float | int]] = {}
    block10_stats: dict[Block, list[float | int]] = {}
    direct = [0.0, 0, 0.0, 0, 0]
    covered_sites = 0
    data_rows = 0
    malformed_rows = 0
    seen_chroms: set[str] = set()
    previous_order: tuple[int, int] | None = None
    previous_coord: tuple[str, int] | None = None
    variants_at_coord: set[tuple[str, str]] = set()
    previous_variant_key: tuple[str, int, str, str] | None = None
    previous_fields: tuple[str, ...] | None = None
    duplicate_rows = 0
    duplicate_stats = [0.0, 0, 0.0, 0, 0]
    source_digest = hashlib.sha256()

    with path.open("rb") as handle:
        header_raw = handle.readline()
        source_digest.update(header_raw)
        header_line = header_raw.decode("utf-8", errors="strict")
        header = tuple(header_line.rstrip("\r\n").split("\t"))
        if header not in EXPECTED_HEADERS:
            raise ValueError(f"Unexpected header in {path}: {header!r}")
        for line_number, raw_line in enumerate(handle, start=2):
            source_digest.update(raw_line)
            line = raw_line.decode("utf-8", errors="strict")
            data_rows += 1
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(header)} fields, "
                    f"observed {len(fields)}"
                )
            try:
                chrom = normalize_chrom(fields[0])
                pos = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
            if chrom not in GRCH37_AUTOSOME_LENGTHS:
                raise ValueError(f"{path}:{line_number}: non-autosomal chromosome {chrom}")
            if not 1 <= pos <= GRCH37_AUTOSOME_LENGTHS[chrom]:
                raise ValueError(f"{path}:{line_number}: position {chrom}:{pos} outside GRCh37")
            order = (int(chrom), pos)
            if previous_order is not None and order < previous_order:
                raise ValueError(
                    f"{path}:{line_number}: CHR/POS order decreased from "
                    f"{previous_order[0]}:{previous_order[1]} to {chrom}:{pos}"
                )
            coord = (chrom, pos)
            if coord != previous_coord:
                variants_at_coord.clear()
                previous_coord = coord
            ref, alt = fields[2].upper(), fields[3].upper()
            variant = (ref, alt)
            variant_key = (chrom, pos, ref, alt)
            row_fields = tuple(fields)
            is_exact_adjacent_duplicate = False
            if variant_key == previous_variant_key:
                if row_fields != previous_fields:
                    raise ValueError(
                        f"{path}:{line_number}: conflicting adjacent duplicate "
                        f"CHR/POS/REF/ALT {chrom}:{pos}:{ref}:{alt}"
                    )
                # Preserve exact adjacent duplicates in primary statistics because
                # the published canonical estimator counted input rows.  Track
                # their contribution separately for a deduplicated diagnostic.
                is_exact_adjacent_duplicate = True
                duplicate_rows += 1
            elif variant in variants_at_coord:
                raise ValueError(
                    f"{path}:{line_number}: non-adjacent duplicate CHR/POS/REF/ALT "
                    f"{chrom}:{pos}:{ref}:{alt}"
                )
            else:
                variants_at_coord.add(variant)
            previous_variant_key = variant_key
            previous_fields = row_fields
            previous_order = order
            seen_chroms.add(chrom)
            try:
                ad_ref = int(fields[5] or 0)
                ad_alt = int(fields[6] or 0)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid allele count") from exc
            if ad_ref < 0 or ad_alt < 0:
                raise ValueError(f"{path}:{line_number}: negative allele count")
            depth = ad_ref + ad_alt
            if depth > 0:
                covered_sites += 1
                direct[4] += 1
                fraction = ad_alt / depth
                if not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
                    raise ValueError(f"{path}:{line_number}: invalid ALT fraction {fraction}")
                add_covered(chrom_stats, chrom)
                add_covered(
                    block5_stats, (chrom, fixed_block_index(chrom, pos, 5_000_000))
                )
                add_covered(
                    block10_stats, (chrom, fixed_block_index(chrom, pos, 10_000_000))
                )
                if is_exact_adjacent_duplicate:
                    duplicate_stats[4] += 1
            if depth < 3:
                continue
            pair = (ref, alt)
            if pair in DAMAGE:
                group = "damage"
            elif (
                ref in "ACGT" and alt in "ACGT" and ref != alt
                and pair not in TRANSITIONS
            ):
                group = "transversion"
            else:
                continue
            # fraction is defined above for every row with nonzero allele depth.
            add_to(chrom_stats, chrom, group, fraction)
            add_to(block5_stats, (chrom, fixed_block_index(chrom, pos, 5_000_000)), group, fraction)
            add_to(block10_stats, (chrom, fixed_block_index(chrom, pos, 10_000_000)), group, fraction)
            if group == "damage":
                direct[0] += fraction
                direct[1] += 1
                if is_exact_adjacent_duplicate:
                    duplicate_stats[0] += fraction
                    duplicate_stats[1] += 1
            else:
                direct[2] += fraction
                direct[3] += 1
                if is_exact_adjacent_duplicate:
                    duplicate_stats[2] += fraction
                    duplicate_stats[3] += 1

    if seen_chroms != set(AUTOSOMES):
        raise ValueError(
            f"{path}: expected exact chromosome set 1..22; "
            f"missing={sorted(set(AUTOSOMES)-seen_chroms, key=int)}, "
            f"extra={sorted(seen_chroms-set(AUTOSOMES))}"
        )
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f"Source changed while being read: {path}")
    return {
        "state": state,
        "source": str(path),
        "size_bytes": before.st_size,
        "mtime_ns": before.st_mtime_ns,
        "source_sha256": source_digest.hexdigest(),
        "data_rows": data_rows,
        "malformed_rows": malformed_rows,
        "exact_adjacent_duplicate_rows": duplicate_rows,
        "duplicate_stats": tuple(duplicate_stats),
        "covered_sites": covered_sites,
        "direct": tuple(direct),
        "chrom": freeze_stats(chrom_stats),
        "5mb": freeze_stats(block5_stats),
        "10mb": freeze_stats(block10_stats),
    }


def add_stats(items: Iterable[Stats]) -> Stats:
    values = list(items)
    return (
        math.fsum(x[0] for x in values),
        sum(x[1] for x in values),
        math.fsum(x[2] for x in values),
        sum(x[3] for x in values),
        sum(x[4] for x in values),
    )


def subtract_stats(total: Stats, removed: Stats) -> Stats:
    result = (
        total[0] - removed[0], total[1] - removed[1],
        total[2] - removed[2], total[3] - removed[3],
        total[4] - removed[4],
    )
    if result[1] < 0 or result[3] < 0 or result[4] < 0:
        raise ArithmeticError("Negative count after block subtraction")
    return result


def estimate(stats: Stats) -> tuple[float, float, float] | None:
    damage_sum, damage_n, tv_sum, tv_n, _covered_n = stats
    if damage_n == 0 or tv_n == 0:
        return None
    damage_pct = 100.0 * damage_sum / damage_n
    tv_pct = 100.0 * tv_sum / tv_n
    return damage_pct, tv_pct, damage_pct - tv_pct


def delta(stats: Stats) -> float:
    result = estimate(stats)
    if result is None:
        raise ValueError("DeltaALT is undefined because a SNP class has zero sites")
    return result[2]


def chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    if chrom in AUTOSOMES:
        return (0, int(chrom))
    return ({"X": (1, 0), "Y": (2, 0), "MT": (3, 0)}).get(chrom, (4, chrom))


def block_sort_key(block: Block) -> tuple[tuple[int, int | str], int]:
    return chrom_sort_key(block[0]), block[1]


def autosomal_blocks(blocks: Mapping[Block, Stats]) -> dict[Block, Stats]:
    return {key: value for key, value in blocks.items() if key[0] in AUTOSOMES}


def state_scope_stats(result: Mapping[str, object], scope: str) -> Stats:
    if scope == "all":
        return result["direct"]  # type: ignore[return-value]
    chrom: Mapping[str, Stats] = result["chrom"]  # type: ignore[assignment]
    if scope == "autosomes":
        return add_stats(chrom.get(c, ZERO_STATS) for c in AUTOSOMES)
    if scope in {"X", "Y"}:
        return chrom.get(scope, ZERO_STATS)
    raise ValueError(scope)


@dataclass(frozen=True)
class Contrast:
    contrast_id: str
    individual: str
    description: str
    terms: tuple[tuple[State, float], ...]


def build_contrasts() -> list[Contrast]:
    contrasts: list[Contrast] = []
    for individual in ("ORT15", "ORT16"):
        def s(library: str, treatment: str) -> State:
            return individual, library, treatment

        contrasts.append(Contrast(
            "raw_nu_minus_raw_fu", individual,
            "Uncorrected non-UDG minus uncorrected full-UDG DeltaALT",
            ((s("nu", "raw"), 1.0), (s("fu", "raw"), -1.0)),
        ))
        for library, treatments in (
            ("fu", ("trim5", "rescale5", "rescale12", "bamrefine5")),
            ("nu", (
                "trim5", "trim10", "rescale5", "rescale10",
                "rescale12", "bamrefine5", "bamrefine10",
            )),
        ):
            for treatment in treatments:
                contrasts.append(Contrast(
                    f"{library}_{treatment}_minus_raw", individual,
                    f"{library} {treatment} minus {library} raw DeltaALT",
                    ((s(library, treatment), 1.0), (s(library, "raw"), -1.0)),
                ))
        contrasts.append(Contrast(
            "nu_bamrefine5_minus_rescale5", individual,
            "Non-UDG bamRefine-5 minus Rescale-5 DeltaALT",
            ((s("nu", "bamrefine5"), 1.0), (s("nu", "rescale5"), -1.0)),
        ))
        for treatment in ("rescale5", "rescale12"):
            contrasts.append(Contrast(
                f"{treatment}_shift_nu_minus_fu", individual,
                f"Difference in {treatment} shift: non-UDG minus full-UDG",
                (
                    (s("nu", treatment), 1.0), (s("nu", "raw"), -1.0),
                    (s("fu", treatment), -1.0), (s("fu", "raw"), 1.0),
                ),
            ))
    assert len(contrasts) == 30
    return contrasts


CONTRASTS = build_contrasts()


def t_critical_975(df: int) -> float:
    # Verified two-sided 95% Student-t critical values for the three fixed
    # block universes. Keeping these release constants avoids a heavy runtime
    # dependency while reproducing the validated Round45 intervals.
    critical = {21: 2.07961384473, 576: 1.96409103161, 288: 1.96823517361}
    try:
        return critical[df]
    except KeyError as exc:
        raise ValueError(f"Unsupported jackknife degrees of freedom: {df}") from exc


def jackknife_summary(full: float, deleted: Sequence[float]) -> dict[str, float | int | bool]:
    n = len(deleted)
    if n < 2:
        raise ValueError("At least two delete-one replicates are required")
    mean_deleted = math.fsum(deleted) / n
    se = math.sqrt((n - 1) / n * math.fsum((x - mean_deleted) ** 2 for x in deleted))
    bias_corrected = n * full - (n - 1) * mean_deleted
    critical = t_critical_975(n - 1)
    lo = bias_corrected - critical * se
    hi = bias_corrected + critical * se
    return {
        "n_blocks": n,
        "full_estimate_pp": full,
        "mean_delete_one_pp": mean_deleted,
        "jackknife_bias_corrected_pp": bias_corrected,
        "jackknife_se_pp": se,
        "t_critical_975": critical,
        "ci95_lo_pp": lo,
        "ci95_hi_pp": hi,
        "ci_excludes_zero": not (lo <= 0.0 <= hi),
    }


def scheme_maps(results: Mapping[State, Mapping[str, object]], scheme: str) -> dict[State, dict[Block, Stats]]:
    output: dict[State, dict[Block, Stats]] = {}
    for state, result in results.items():
        if scheme == "loco":
            chrom: Mapping[str, Stats] = result["chrom"]  # type: ignore[assignment]
            missing = [c for c in AUTOSOMES if c not in chrom]
            if missing:
                raise ValueError(f"{state_label(state)} lacks eligible sites on chromosomes {missing}")
            output[state] = {(c, 0): chrom[c] for c in AUTOSOMES}
        else:
            observed = autosomal_blocks(result[scheme])  # type: ignore[arg-type]
            universe = COMMON_5MB_BLOCKS if scheme == "5mb" else COMMON_10MB_BLOCKS
            unexpected = set(observed) - set(universe)
            if unexpected:
                raise ValueError(f"Unexpected {scheme} blocks for {state_label(state)}: {unexpected}")
            output[state] = {block: observed.get(block, ZERO_STATS) for block in universe}
    return output


def compute_jackknives(results: Mapping[State, Mapping[str, object]]) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    state_summaries: list[dict] = []
    state_replicates: list[dict] = []
    contrast_summaries: list[dict] = []
    contrast_replicates: list[dict] = []
    for scheme in ("loco", "5mb", "10mb"):
        maps = scheme_maps(results, scheme)
        for state in STATES:
            blocks = maps[state]
            ordered = sorted(blocks, key=block_sort_key)
            full_stats = add_stats(blocks[b] for b in ordered)
            full = delta(full_stats)
            full_components = estimate(full_stats)
            assert full_components is not None
            deleted: list[float] = []
            for block in ordered:
                value = delta(subtract_stats(full_stats, blocks[block]))
                deleted.append(value)
                state_replicates.append({
                    "scheme": scheme, "individual": state[0], "library": state[1],
                    "treatment": state[2], "omitted_chrom": block[0],
                    "omitted_block_index_0based": block[1], "delete_one_deltaalt_pp": value,
                })
            state_summaries.append({
                "scheme": scheme, "individual": state[0], "library": state[1],
                "treatment": state[2],
                "damage_alt_fraction_sum": full_stats[0], "damage_n_dp3": full_stats[1],
                "transversion_alt_fraction_sum": full_stats[2],
                "transversion_n_dp3": full_stats[3], "covered_sites": full_stats[4],
                "damage_alt_pct": full_components[0],
                "transversion_alt_pct": full_components[1],
                "estimate_unit": "percentage_points",
                **jackknife_summary(full, deleted),
            })

        for contrast in CONTRASTS:
            # Every state map uses the same GRCh37 block universe, so B and the
            # omitted physical region are identical for all states and contrasts.
            universe = sorted(maps[contrast.terms[0][0]], key=block_sort_key)
            if any(set(maps[state]) != set(universe) for state, _ in contrast.terms):
                raise AssertionError("Paired contrast states do not share a block universe")
            totals = {
                state: add_stats(maps[state].get(block, ZERO_STATS) for block in universe)
                for state, _ in contrast.terms
            }
            full = math.fsum(coef * delta(totals[state]) for state, coef in contrast.terms)
            deleted_values: list[float] = []
            for block in universe:
                value = math.fsum(
                    coef * delta(subtract_stats(totals[state], maps[state].get(block, ZERO_STATS)))
                    for state, coef in contrast.terms
                )
                deleted_values.append(value)
                contrast_replicates.append({
                    "scheme": scheme, "individual": contrast.individual,
                    "contrast_id": contrast.contrast_id, "omitted_chrom": block[0],
                    "omitted_block_index_0based": block[1], "delete_one_contrast_pp": value,
                })
            contrast_summaries.append({
                "scheme": scheme, "individual": contrast.individual,
                "contrast_id": contrast.contrast_id, "description": contrast.description,
                "terms_json": json.dumps([
                    {"state": state_label(state), "coefficient": coef}
                    for state, coef in contrast.terms
                ], separators=(",", ":")),
                "estimate_unit": "percentage_points",
                **jackknife_summary(full, deleted_values),
            })
    return state_summaries, state_replicates, contrast_summaries, contrast_replicates


def read_canonical(path: Path) -> dict[State, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    data = {(r["individual"], r["library"], r["treatment"]): r for r in rows}
    if len(rows) != 26 or len(data) != 26 or set(data) != set(STATES):
        raise ValueError("Canonical CSV must contain exactly the expected 26 unique states")
    return data


def validate_canonical(results: Mapping[State, Mapping[str, object]], canonical_path: Path) -> tuple[list[dict], bool]:
    canonical = read_canonical(canonical_path)
    rows: list[dict] = []
    all_pass = True
    for state in STATES:
        result = results[state]
        observed = estimate(result["direct"])  # type: ignore[arg-type]
        if observed is None:
            raise ValueError(f"Undefined full-scope DeltaALT for {state_label(state)}")
        ts, tv, diff = observed
        expected = canonical[state]
        float_checks = {
            "ts_alt_pct": (f"{ts:.6f}", expected["ts_alt_pct"]),
            "tv_alt_pct": (f"{tv:.6f}", expected["tv_alt_pct"]),
            "ts_minus_tv_pp": (f"{diff:+.6f}", expected["ts_minus_tv_pp"]),
        }
        direct: Stats = result["direct"]  # type: ignore[assignment]
        int_checks = {
            "ts_n_dp3": (direct[1], int(expected["ts_n_dp3"])),
            "tv_n_dp3": (direct[3], int(expected["tv_n_dp3"])),
            "covered_sites": (int(result["covered_sites"]), int(expected["covered_sites"])),
        }
        raw_covered = int(results[(state[0], state[1], "raw")]["covered_sites"])
        retention = 100.0 * int(result["covered_sites"]) / raw_covered
        float_checks["retention_pct"] = (f"{retention:.6f}", expected["retention_pct"])
        source_ok = Path(str(result["source"])).name == expected["source_file"]
        malformed_ok = int(result["malformed_rows"]) == 0
        passed = (
            all(obs == exp for obs, exp in float_checks.values())
            and all(obs == exp for obs, exp in int_checks.values())
            and source_ok and malformed_ok
        )
        all_pass = all_pass and passed
        rows.append({
            "individual": state[0], "library": state[1], "treatment": state[2],
            "computed_ts_alt_pct": float_checks["ts_alt_pct"][0],
            "canonical_ts_alt_pct": float_checks["ts_alt_pct"][1],
            "computed_tv_alt_pct": float_checks["tv_alt_pct"][0],
            "canonical_tv_alt_pct": float_checks["tv_alt_pct"][1],
            "computed_deltaalt_pp": float_checks["ts_minus_tv_pp"][0],
            "canonical_deltaalt_pp": float_checks["ts_minus_tv_pp"][1],
            "computed_ts_n_dp3": direct[1], "canonical_ts_n_dp3": int_checks["ts_n_dp3"][1],
            "computed_tv_n_dp3": direct[3], "canonical_tv_n_dp3": int_checks["tv_n_dp3"][1],
            "computed_covered_sites": result["covered_sites"],
            "canonical_covered_sites": int_checks["covered_sites"][1],
            "computed_retention_pct": float_checks["retention_pct"][0],
            "canonical_retention_pct": float_checks["retention_pct"][1],
            "source_basename_match": source_ok, "malformed_rows": result["malformed_rows"],
            "exact_adjacent_duplicate_rows_preserved": result["exact_adjacent_duplicate_rows"],
            "duplicate_covered_sites": result["duplicate_stats"][4],
            "duplicate_damage_n_dp3": result["duplicate_stats"][1],
            "duplicate_transversion_n_dp3": result["duplicate_stats"][3],
            "validation_pass": passed,
        })
    return rows, all_pass


def check_internal_consistency(results: Mapping[State, Mapping[str, object]]) -> None:
    for state in STATES:
        result = results[state]
        chrom: Mapping[str, Stats] = result["chrom"]  # type: ignore[assignment]
        by_chrom = add_stats(chrom.values())
        direct: Stats = result["direct"]  # type: ignore[assignment]
        if (by_chrom[1], by_chrom[3], by_chrom[4]) != (direct[1], direct[3], direct[4]):
            raise AssertionError(f"Chromosome counts do not sum to direct totals for {state_label(state)}")
        sum_tolerance = max(1e-8, 2e-12 * max(abs(direct[0]), abs(direct[2]), 1.0))
        if not math.isclose(by_chrom[0], direct[0], rel_tol=0.0, abs_tol=sum_tolerance):
            raise AssertionError(f"Damage sums disagree for {state_label(state)}")
        if not math.isclose(by_chrom[2], direct[2], rel_tol=0.0, abs_tol=sum_tolerance):
            raise AssertionError(f"Transversion sums disagree for {state_label(state)}")
        for scheme in ("5mb", "10mb"):
            combined = add_stats(result[scheme].values())  # type: ignore[union-attr]
            if (combined[1], combined[3], combined[4]) != (direct[1], direct[3], direct[4]):
                raise AssertionError(f"{scheme} counts do not sum for {state_label(state)}")
            if not math.isclose(combined[0], direct[0], rel_tol=0.0, abs_tol=sum_tolerance):
                raise AssertionError(f"{scheme} damage sums disagree for {state_label(state)}")
            if not math.isclose(combined[2], direct[2], rel_tol=0.0, abs_tol=sum_tolerance):
                raise AssertionError(f"{scheme} transversion sums disagree for {state_label(state)}")


def fmt(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.12g}"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    if fieldnames is None:
        if not rows:
            raise ValueError(f"Cannot infer columns for empty output {path.name}")
        fieldnames = list(rows[0])
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fieldnames})
    os.replace(tmp, path)


def sufficient_rows(results: Mapping[State, Mapping[str, object]], scheme: str) -> list[dict]:
    rows: list[dict] = []
    for state in STATES:
        blocks: Mapping[object, Stats] = results[state][scheme]  # type: ignore[assignment]
        if scheme == "chrom":
            ordered = list(AUTOSOMES)
            iterator = ((chrom, 0, blocks[chrom]) for chrom in ordered)
            size = None
        else:
            size = 5_000_000 if scheme == "5mb" else 10_000_000
            universe = COMMON_5MB_BLOCKS if scheme == "5mb" else COMMON_10MB_BLOCKS
            iterator = (
                (block[0], block[1], blocks.get(block, ZERO_STATS)) for block in universe
            )  # type: ignore[union-attr]
        for chrom, block_index, stats in iterator:
            est = estimate(stats)
            if size is None:
                block_start, block_end = "", ""
            else:
                block_start, block_end = fixed_block_bounds(chrom, block_index, size)
            rows.append({
                "individual": state[0], "library": state[1], "treatment": state[2],
                "scheme": scheme, "chrom": chrom,
                "block_id": (
                    f"chr{chrom}" if size is None
                    else f"chr{chrom}:{block_start}-{block_end}"
                ),
                "block_index_0based": "" if size is None else block_index,
                "block_start_1based": block_start,
                "block_end_inclusive": block_end,
                "damage_alt_fraction_sum": stats[0], "damage_n_dp3": stats[1],
                "transversion_alt_fraction_sum": stats[2], "transversion_n_dp3": stats[3],
                "covered_sites": stats[4],
                "damage_alt_pct": "" if est is None else est[0],
                "transversion_alt_pct": "" if est is None else est[1],
                "deltaalt_pp": "" if est is None else est[2],
                "deltaalt_unit": "percentage_points",
            })
    return rows


def scope_rows(results: Mapping[State, Mapping[str, object]]) -> list[dict]:
    rows: list[dict] = []
    for state in STATES:
        for scope in ("all", "autosomes", "X", "Y"):
            stats = state_scope_stats(results[state], scope)
            est = estimate(stats)
            rows.append({
                "individual": state[0], "library": state[1], "treatment": state[2],
                "scope": scope, "damage_alt_fraction_sum": stats[0],
                "damage_n_dp3": stats[1], "transversion_alt_fraction_sum": stats[2],
                "transversion_n_dp3": stats[3], "covered_sites": stats[4],
                "damage_alt_pct": "" if est is None else est[0],
                "transversion_alt_pct": "" if est is None else est[1],
                "deltaalt_pp": "" if est is None else est[2],
                "deltaalt_unit": "percentage_points",
                "estimate_defined": est is not None,
            })
    return rows


def duplicate_diagnostic_rows(results: Mapping[State, Mapping[str, object]]) -> list[dict]:
    rows: list[dict] = []
    for state in STATES:
        primary: Stats = results[state]["direct"]  # type: ignore[assignment]
        duplicate: Stats = results[state]["duplicate_stats"]  # type: ignore[assignment]
        deduplicated = subtract_stats(primary, duplicate)
        primary_est = estimate(primary)
        dedup_est = estimate(deduplicated)
        assert primary_est is not None and dedup_est is not None
        rows.append({
            "individual": state[0], "library": state[1], "treatment": state[2],
            "exact_adjacent_duplicate_rows": results[state]["exact_adjacent_duplicate_rows"],
            "duplicate_covered_sites": duplicate[4],
            "duplicate_damage_alt_fraction_sum": duplicate[0],
            "duplicate_damage_n_dp3": duplicate[1],
            "duplicate_transversion_alt_fraction_sum": duplicate[2],
            "duplicate_transversion_n_dp3": duplicate[3],
            "published_row_counting_deltaalt_pp": primary_est[2],
            "deduplicated_diagnostic_deltaalt_pp": dedup_est[2],
            "deduplicated_minus_published_pp": dedup_est[2] - primary_est[2],
            "published_covered_sites": primary[4],
            "deduplicated_covered_sites": deduplicated[4],
            "diagnostic_only": "yes; not used in primary jackknife",
        })
    return rows


def run_self_tests() -> None:
    """Fast deterministic design checks; source-content gates run during streaming."""
    probe_root = Path("/external/deltaalt-source-root")
    assert len(STATES) == 26 and len({source_path(probe_root, s) for s in STATES}) == 26
    assert len(CONTRASTS) == 30
    for size, universe in (
        (5_000_000, COMMON_5MB_BLOCKS), (10_000_000, COMMON_10MB_BLOCKS),
    ):
        assert len(universe) == len(set(universe))
        for chrom in AUTOSOMES:
            blocks = [block for block in universe if block[0] == chrom]
            bounds = [fixed_block_bounds(block[0], block[1], size) for block in blocks]
            assert bounds[0][0] == 1
            assert bounds[-1][1] == GRCH37_AUTOSOME_LENGTHS[chrom]
            assert all(left[1] + 1 == right[0] for left, right in zip(bounds, bounds[1:]))
            assert all(end - start + 1 >= size / 2 for start, end in bounds)
    probe: Stats = (3.0, 2, 1.0, 2, 9)
    assert delta(probe) == 100.0
    assert subtract_stats(probe, ZERO_STATS) == probe
    try:
        subtract_stats(ZERO_STATS, probe)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("Negative-count subtraction control did not fail")


def write_validation_report(path: Path, results: Mapping[State, Mapping[str, object]]) -> None:
    x_present = any("X" in result["chrom"] for result in results.values())  # type: ignore[operator]
    y_present = any("Y" in result["chrom"] for result in results.values())  # type: ignore[operator]
    duplicate_total = sum(int(r["exact_adjacent_duplicate_rows"]) for r in results.values())
    lines = [
        "# DeltaALT uncertainty validation\n",
        "- Canonical 26-state validation: **PASS** (exact six-decimal strings and integer counts).",
        "- Source access: read-only streaming; SHA-256 computed from raw bytes in the same pass.",
        "- Row gates: exact nine-field schema, nonnegative counts, finite recomputed fractions,",
        "  monotonic coordinates, and conflicting/non-adjacent duplicate rejection.",
        f"- Exact-identical adjacent duplicate rows preserved for canonical reproduction: {duplicate_total}.",
        "  Their contributions and a deduplicated diagnostic are reported separately.",
        "- Chromosome gate: exact set 1–22 in every source TSV.",
        f"- Five-megabase physical block universe: {len(COMMON_5MB_BLOCKS)} GRCh37 blocks.",
        f"- Ten-megabase physical block universe: {len(COMMON_10MB_BLOCKS)} GRCh37 blocks.",
        "- Short terminal remainder: merged into the preceding block when shorter than half",
        "  the target width.",
        "- Paired deletion: the same full-universe physical block is omitted from every",
        "  term in a contrast.",
        "- Site sets: correction-specific, as prespecified; no common-callable intersection.",
        "- Common-callable sensitivity analysis: **not performed by this script**; it would",
        "  require a separate synchronized multi-file scan and is outside this run's scope.",
        f"- X present: {'yes' if x_present else 'no'}; Y present: {'yes' if y_present else 'no'}.",
        "  X/Y estimates are therefore reported as undefined when those chromosomes are absent.",
        "- Scientific source files modified: **no**.\n",
    ]
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    os.replace(tmp, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def output_directory(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="New or empty output directory")
    parser.add_argument(
        "--source-root", type=Path, required=True,
        help="External GLIMPSE analysis root containing the 26 per-site TSVs",
    )
    parser.add_argument("--canonical-csv", type=Path, default=DEFAULT_CANONICAL)
    parser.add_argument("--workers", type=int, default=8, help="Parallel readers (1-8; default 8)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_self_tests()
    if not 1 <= args.workers <= 8:
        raise ValueError("--workers must be between 1 and 8")
    if not args.canonical_csv.is_file():
        raise FileNotFoundError(f"Canonical CSV not found: {args.canonical_csv}")
    tasks = [(state, str(source_path(args.source_root, state))) for state in STATES]
    missing = [path for _, path in tasks if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError("Missing source TSV(s):\n" + "\n".join(missing))
    output_directory(args.output_dir)
    started = utc_now()
    results: dict[State, Mapping[str, object]] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(stream_one, task): task[0] for task in tasks}
        for future in as_completed(futures):
            state = futures[future]
            results[state] = future.result()
            print(f"completed {state_label(state)}", file=sys.stderr, flush=True)
    if set(results) != set(STATES):
        raise AssertionError("Did not receive all 26 worker results")

    check_internal_consistency(results)
    validation, validation_pass = validate_canonical(results, args.canonical_csv)
    write_csv(args.output_dir / "canonical_validation.csv", validation)
    if not validation_pass:
        raise RuntimeError(
            f"Canonical validation failed; inspect {args.output_dir / 'canonical_validation.csv'}"
        )

    write_csv(args.output_dir / "sufficient_stats_chromosome.csv", sufficient_rows(results, "chrom"))
    write_csv(args.output_dir / "sufficient_stats_5mb.csv", sufficient_rows(results, "5mb"))
    write_csv(args.output_dir / "sufficient_stats_10mb.csv", sufficient_rows(results, "10mb"))
    write_csv(args.output_dir / "state_estimates_by_scope.csv", scope_rows(results))
    write_csv(
        args.output_dir / "exact_duplicate_diagnostic.csv",
        duplicate_diagnostic_rows(results),
    )

    state_summary, state_rep, contrast_summary, contrast_rep = compute_jackknives(results)
    write_csv(args.output_dir / "jackknife_state_summary.csv", state_summary)
    write_csv(args.output_dir / "jackknife_state_delete_one.csv", state_rep)
    write_csv(args.output_dir / "jackknife_contrast_summary.csv", contrast_summary)
    write_csv(args.output_dir / "jackknife_contrast_delete_one.csv", contrast_rep)

    source_rows = []
    for state in STATES:
        result = results[state]
        source_rows.append({
            "individual": state[0], "library": state[1], "treatment": state[2],
            "source_file": Path(str(result["source"])).name,
            "size_bytes": result["size_bytes"],
            "source_sha256": result["source_sha256"],
            "data_rows": result["data_rows"],
            "covered_sites": result["covered_sites"], "malformed_rows": result["malformed_rows"],
            "exact_adjacent_duplicate_rows": result["exact_adjacent_duplicate_rows"],
            "duplicate_covered_sites": result["duplicate_stats"][4],
            "duplicate_damage_n_dp3": result["duplicate_stats"][1],
            "duplicate_transversion_n_dp3": result["duplicate_stats"][3],
        })
    write_csv(args.output_dir / "source_manifest.csv", source_rows)
    write_validation_report(args.output_dir / "VALIDATION.md", results)

    script = Path(__file__).resolve()
    manifest = {
        "status": "complete",
        "started_utc": started,
        "finished_utc": utc_now(),
        "platform": platform.platform(),
        "python": sys.version,
        "script": script.name,
        "script_sha256": sha256(script),
        "canonical_csv": args.canonical_csv.name,
        "canonical_csv_sha256": sha256(args.canonical_csv),
        "source_root": "external; paths intentionally not serialized",
        "source_sha256_by_state": {
            state_label(state): results[state]["source_sha256"] for state in STATES
        },
        "workers": args.workers,
        "states": len(STATES),
        "block_schemes": [
            "leave-one-autosome-out",
            "5mb-autosomal-GRCh37-short-terminal-remainder-merged-if-under-2.5mb",
            "10mb-autosomal-GRCh37-short-terminal-remainder-merged-if-under-5mb",
        ],
        "site_sets": "correction-specific; no intersection between states",
        "paired_deletion": "same genomic block omitted from every term of a contrast",
        "deltaalt_definition": (
            "100 * (mean[AD_ALT/(AD_REF+AD_ALT)] at C>T/G>A SNPs - "
            "mean at transversion SNPs), with AD_REF+AD_ALT >= 3"
        ),
        "canonical_validation_pass": True,
        "internal_self_tests_pass": True,
        "fixed_block_counts": {
            "5mb": len(COMMON_5MB_BLOCKS), "10mb": len(COMMON_10MB_BLOCKS)
        },
        "common_callable_sensitivity": "not performed; outside this correction-specific run",
        "exact_adjacent_duplicate_rows_preserved": sum(
            int(results[state]["exact_adjacent_duplicate_rows"]) for state in STATES
        ),
        "duplicate_policy": (
            "preserve exact-identical adjacent duplicate input rows in canonical/primary "
            "statistics; reject conflicting or non-adjacent duplicate variant keys; "
            "emit a separate deduplicated diagnostic"
        ),
        "source_files_modified": False,
    }
    tmp = args.output_dir / "run_manifest.json.tmp"
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, args.output_dir / "run_manifest.json")
    print(f"complete: {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
