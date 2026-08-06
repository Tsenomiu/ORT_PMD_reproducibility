#!/usr/bin/env python3
"""Exact conservative allele harmonization for reference, imputed, and PHAP BIMs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
MISSING = {"0", "."}


def read_bim(path: Path):
    result = {}
    with path.open() as handle:
        for raw in handle:
            chrom, identifier, _, pos, a1, a2 = raw.split()[:6]
            result[identifier] = (chrom, pos, a1.upper(), a2.upper())
    return result


def classify(alleles: set[str], reference: set[str]):
    observed = alleles - MISSING; incomplete = bool(alleles & MISSING) or len(observed) < 2
    if not observed: return "missing", None, incomplete
    if observed <= reference: return "same", False, incomplete
    if {COMP.get(x, x) for x in observed} <= reference: return "flip", True, incomplete
    return "mismatch", None, incomplete


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path); parser.add_argument("imputed", type=Path)
    parser.add_argument("pseudohaploid", type=Path); parser.add_argument("output_dir", type=Path)
    args = parser.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    ref, imp, phap = map(read_bim, (args.reference, args.imputed, args.pseudohaploid))
    common = sorted(set(ref) & set(imp) & set(phap), key=lambda x: tuple(map(int, x.split(":"))))
    kept = []; flip_imp = []; flip_phap = []; stats = {"n_common_raw": len(common), "n_ambiguous": 0, "n_mismatch": 0, "n_query_monomorphic": 0}
    for identifier in common:
        _, _, r1, r2 = ref[identifier]; _, _, i1, i2 = imp[identifier]; _, _, p1, p2 = phap[identifier]
        if (r1, r2) in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}:
            stats["n_ambiguous"] += 1; continue
        ci, cq = classify({i1, i2}, {r1, r2}), classify({p1, p2}, {r1, r2})
        if ci[0] == "mismatch" or cq[0] == "mismatch": stats["n_mismatch"] += 1; continue
        if ci[2] or cq[2] or ci[0] == "missing" or cq[0] == "missing": stats["n_query_monomorphic"] += 1; continue
        kept.append(identifier)
        if ci[1]: flip_imp.append(identifier)
        if cq[1]: flip_phap.append(identifier)
    (args.output_dir / "common_sites_matched48.txt").write_text("\n".join(kept) + "\n")
    (args.output_dir / "flip_imp.txt").write_text("\n".join(flip_imp) + ("\n" if flip_imp else ""))
    (args.output_dir / "flip_phap.txt").write_text("\n".join(flip_phap) + ("\n" if flip_phap else ""))
    stats.update({"n_kept": len(kept), "n_flip_imp": len(flip_imp), "n_flip_phap": len(flip_phap),
                  "n_dropped_total": stats["n_ambiguous"] + stats["n_mismatch"] + stats["n_query_monomorphic"]})
    (args.output_dir / "site_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

