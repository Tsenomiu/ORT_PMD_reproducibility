#!/usr/bin/env python3
"""Add/replace FORMAT/DS using DS=P(0/1)+2*P(1/1), rounded as executed."""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=Path)
    args = parser.parse_args()
    inserted = False
    with open_text(args.vcf) as handle:
        for raw in handle:
            if raw.startswith("##FORMAT=<ID=DS,"):
                inserted = True
            if raw.startswith("#CHROM") and not inserted:
                print('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage reconstructed as GP[1]+2*GP[2]">')
            if raw.startswith("#"):
                sys.stdout.write(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            keys = fields[8].split(":")
            if "GP" not in keys:
                raise ValueError(f"FORMAT/GP absent at {fields[0]}:{fields[1]}")
            gp_i = keys.index("GP")
            if "DS" in keys:
                ds_i = keys.index("DS")
            else:
                keys.append("DS"); ds_i = len(keys) - 1
            for i in range(9, len(fields)):
                values = fields[i].split(":")
                values += ["."] * (len(keys) - len(values))
                gp = values[gp_i].split(",")
                if len(gp) != 3 or "." in gp:
                    values[ds_i] = "."
                else:
                    values[ds_i] = f"{float(gp[1]) + 2.0 * float(gp[2]):.4f}"
                fields[i] = ":".join(values)
            fields[8] = ":".join(keys)
            print("\t".join(fields))


if __name__ == "__main__":
    main()
