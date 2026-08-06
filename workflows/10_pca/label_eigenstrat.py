#!/usr/bin/env python3
"""Label reference populations and mark IMP__/PHAP__ samples as projected."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ind", type=Path); parser.add_argument("iid_population", type=Path); parser.add_argument("output", type=Path)
    args = parser.parse_args()
    mapping = {}
    with args.iid_population.open() as handle:
        for raw in handle:
            fields = raw.rstrip().split("\t")
            if len(fields) >= 2: mapping[fields[0]] = fields[1]
    output = []; n_ref = n_query = 0
    with args.ind.open() as handle:
        for raw in handle:
            fields = raw.split(); iid = fields[0]; sex = fields[1] if len(fields) > 1 else "U"
            if iid.startswith("IMP__"): group = "Projected_IMP"; n_query += 1
            elif iid.startswith("PHAP__"): group = "Projected_PHAP"; n_query += 1
            elif iid in mapping: group = mapping[iid]; n_ref += 1
            else: raise ValueError(f"No reference population for {iid}")
            output.append(f"{iid}\t{sex}\t{group}\n")
    if n_ref != 985 or n_query != 48: raise ValueError(f"Expected 985 reference + 48 query; saw {n_ref} + {n_query}")
    args.output.write_text("".join(output))


if __name__ == "__main__":
    main()

