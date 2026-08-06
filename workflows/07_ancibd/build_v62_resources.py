#!/usr/bin/env python3
"""Build the exact AADR-v62/1000G-AF/common-VCF ancIBD resource intersection."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

AUTOSOMES = range(1, 23)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v62", required=True, type=Path)
    parser.add_argument("--af-pattern", required=True)
    parser.add_argument("--vcf-sites-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-markers", type=int, default=1_100_313)
    parser.add_argument("--expected-span-cm", type=float, default=3538.8511)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in ("markers", "afs", "map_by_chromosome", "map"):
        (args.output / name).mkdir(exist_ok=True)

    v62: dict[int, dict[int, tuple[str, float, str, str, str]]] = {c: {} for c in AUTOSOMES}
    map_lines = []
    with args.v62.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.split()
            if len(fields) < 6:
                raise ValueError(f"Malformed AADR row {line_number}")
            snp, chrom_text, map_text, pos_text, ref, alt = fields[:6]
            try:
                chrom = int(chrom_text)
            except ValueError:
                continue
            if chrom not in v62:
                continue
            pos, genetic_map = int(pos_text), float(map_text)
            if pos in v62[chrom]:
                raise ValueError(f"Duplicate AADR position chr{chrom}:{pos}")
            v62[chrom][pos] = (snp, genetic_map, ref, alt, raw if raw.endswith("\n") else raw + "\n")
            map_lines.append(v62[chrom][pos][4])
    (args.output / "map/v62.autosomes.Morgan.snp").write_text("".join(map_lines))

    manifest = {"schema": "ancIBD-v62-exact-intersection-v1", "join_key": ["chromosome", "position", "REF", "ALT"], "chromosomes": {}}
    total = 0; span_m = 0.0
    for chrom in AUTOSOMES:
        catalog_path = args.vcf_sites_dir / f"vcf_sites_ch{chrom}.tsv"
        catalog = {}
        with catalog_path.open() as handle:
            for raw in handle:
                c, pos, ref, alt = raw.split()[:4]
                if int(c) != chrom or int(pos) in catalog:
                    raise ValueError(f"Invalid catalog row: {raw.rstrip()}")
                catalog[int(pos)] = (ref, alt)
        af_path = Path(args.af_pattern.replace("{CHROM}", str(chrom)))
        joined = []
        with af_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"chr", "pos", "ref", "alt", "af"}.issubset(reader.fieldnames or []):
                raise ValueError(f"Unexpected AF header: {af_path}")
            for row in reader:
                pos = int(row["pos"]); aadr = v62[chrom].get(pos)
                if aadr and (row["ref"], row["alt"]) == aadr[2:4] == catalog.get(pos):
                    frequency = float(row["af"])
                    if not 0 <= frequency <= 1:
                        raise ValueError(f"Invalid AF chr{chrom}:{pos}")
                    joined.append((pos, row["ref"], row["alt"], row["af"], aadr[1]))
        joined.sort()
        if len(joined) < 10_000:
            raise ValueError(f"Implausibly few joined sites on chromosome {chrom}")
        marker = args.output / f"markers/v62_targets_ch{chrom}.tsv"
        af_out = args.output / f"afs/v62_1000G_AF_ch{chrom}.tsv"
        map_out = args.output / f"map_by_chromosome/v62_map_ch{chrom}.tsv"
        marker.write_text("".join(f"{chrom}\t{x[0]}\n" for x in joined))
        af_out.write_text("chr\tpos\tref\talt\taf\n" + "".join(f"{chrom}\t{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\n" for x in joined))
        map_out.write_text("".join(f"{x[0]}\t{x[4]:.8f}\n" for x in joined))
        first_m, last_m = joined[0][4], joined[-1][4]
        span_m += last_m - first_m; total += len(joined)
        manifest["chromosomes"][str(chrom)] = {"markers": len(joined), "canonical_first_M": first_m, "canonical_last_M": last_m}
    span_cm = 100 * span_m
    manifest["exact_join_markers"] = total
    manifest["canonical_callable_first_to_last_span_cM"] = span_cm
    if total != args.expected_markers:
        raise ValueError(f"Expected {args.expected_markers} markers, observed {total}")
    if abs(span_cm - args.expected_span_cm) > 0.0002:
        raise ValueError(f"Expected {args.expected_span_cm:.4f} cM, observed {span_cm:.4f}")
    (args.output / "RESOURCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"markers": total, "span_cM": span_cm}))


if __name__ == "__main__":
    main()

