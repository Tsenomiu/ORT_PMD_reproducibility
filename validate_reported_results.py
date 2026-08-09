#!/usr/bin/env python3
"""Validate compact kinship and mitochondrial tables reported in the study."""

from __future__ import annotations

import csv
import math
import statistics
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def close(observed: float, expected: float, tolerance: float = 1e-12) -> None:
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"{observed} != {expected}")


def validate_tkgwv2() -> None:
    data = rows(ROOT / "data/summary/kinship/tkgwv2_pair_results.tsv")
    if len(data) != 16:
        raise AssertionError(f"expected 16 TKGWV2 states, found {len(data)}")
    by_pair = {(r["sample1"], r["sample2"]): r for r in data}
    ort_pair = by_pair[("ORT15_fu_bamrefine5", "ORT16_fu_bamrefine5")]
    if int(ort_pair["used_snps"]) != 4517214:
        raise AssertionError("unexpected plotted TKGWV2 SNP count")
    close(float(ort_pair["HRC"]), 0.2356)
    if Counter(r["classification"] for r in data) != {
        "first_degree": 11,
        "second_degree": 5,
    }:
        raise AssertionError("unexpected TKGWV2 classification counts")


def validate_readv2() -> None:
    data = rows(ROOT / "data/summary/kinship/readv2_cemetery_210_pairs.tsv")
    if len(data) != 210:
        raise AssertionError(f"expected 210 READv2 pairs, found {len(data)}")
    expected = {
        "Unrelated": 198,
        "First Degree": 3,
        "Second Degree": 4,
        "Third Degree": 5,
    }
    if Counter(r["Rel"] for r in data) != expected:
        raise AssertionError("unexpected READv2 relationship counts")
    ort_pair = next(r for r in data if r["PairIndividuals"] == "ORT15,ORT16")
    close(float(ort_pair["KinshipCoefficient"]), 0.23096666189772574)
    if int(ort_pair["OverlapNSNPs"]) != 169172 or ort_pair["1st_Type"] != "Parent-offspring":
        raise AssertionError("unexpected ORT15--ORT16 READv2 result")
    related = [r for r in data if r["Rel"] != "Unrelated"]
    if min(int(r["OverlapNSNPs"]) for r in related) != 143526:
        raise AssertionError("unexpected minimum related-pair overlap")
    second = [float(r["KinshipCoefficient"]) for r in data if r["Rel"] == "Second Degree"]
    third = [float(r["KinshipCoefficient"]) for r in data if r["Rel"] == "Third Degree"]
    if (round(min(second), 3), round(max(second), 3)) != (0.097, 0.112):
        raise AssertionError("second-degree range changed")
    if (round(min(third), 3), round(max(third), 3)) != (0.048, 0.090):
        raise AssertionError("third-degree range changed")
    low = [
        float(r["KinshipCoefficient"])
        for r in data
        if r["Rel"] == "Unrelated" and int(r["OverlapNSNPs"]) < 200000
    ]
    high = [
        float(r["KinshipCoefficient"])
        for r in data
        if r["Rel"] == "Unrelated" and int(r["OverlapNSNPs"]) >= 200000
    ]
    if round(statistics.median(low), 3) != -0.003:
        raise AssertionError("low-overlap unrelated median changed")
    if round(statistics.median(high), 3) != 0.001:
        raise AssertionError("high-overlap unrelated median changed")

    reported = rows(ROOT / "data/summary/kinship/readv2_ort15_ort16.tsv")
    if {r["analysis"] for r in reported} != {"primary_chr1_22_X_Y", "autosome_only_sensitivity"}:
        raise AssertionError("unexpected READv2 reported-analysis rows")
    primary = next(r for r in reported if r["analysis"] == "primary_chr1_22_X_Y")
    for field, cohort_field in (
        ("PairIndividuals", "PairIndividuals"),
        ("Rel", "Rel"),
        ("OverlapNSNPs", "OverlapNSNPs"),
        ("KinshipCoefficient", "KinshipCoefficient"),
        ("First_degree_subtype", "1st_Type"),
    ):
        if primary[field] != ort_pair[cohort_field]:
            raise AssertionError(f"READv2 primary view differs from cohort table in {field}")
    sensitivity = next(r for r in reported if r["analysis"] == "autosome_only_sensitivity")
    if int(sensitivity["OverlapNSNPs"]) != 157451:
        raise AssertionError("unexpected autosome-only READv2 overlap")
    close(float(sensitivity["KinshipCoefficient"]), 0.23159213948114488)
    if sensitivity["First_degree_subtype"] != "Parent-offspring":
        raise AssertionError("unexpected autosome-only READv2 subtype")


def validate_king_ibs0() -> None:
    data = {r["analysis"]: r for r in rows(ROOT / "data/summary/kinship/king_ibs0_summary.tsv")}
    if set(data) != {"genomewide", "chromosome1_high_confidence"}:
        raise AssertionError("unexpected KING/IBS0 summary rows")
    genomewide = data["genomewide"]
    if int(genomewide["sites"]) != 5272558 or int(genomewide["IBS0"]) != 17988:
        raise AssertionError("unexpected genome-wide KING/IBS0 counts")
    close(float(genomewide["KING_robust"]), 0.2379)
    close(float(genomewide["observed_IBS0_per_site"]), 0.00341)
    close(float(genomewide["unrelated_expected"]), 0.07341)
    close(float(genomewide["full_sibling_expected"]), 0.01835)
    high_confidence = data["chromosome1_high_confidence"]
    if int(high_confidence["sites"]) != 324285 or int(high_confidence["IBS0"]) != 36:
        raise AssertionError("unexpected high-confidence chromosome-1 IBS0 counts")
    close(float(high_confidence["observed_IBS0_per_site"]), 0.00011)


def validate_mtdna() -> None:
    calls = rows(ROOT / "data/summary/mtdna/haplogrep_classification.tsv")
    if len(calls) != 2 or {r["haplogroup"] for r in calls} != {"D4o1"}:
        raise AssertionError("unexpected mitochondrial classifications")
    scores = {r["sample"]: float(r["quality_score"]) for r in calls}
    close(scores["ORT15"], 0.8911)
    close(scores["ORT16"], 0.9107)
    concordance = rows(ROOT / "data/summary/mtdna/mtdna_pair_concordance.tsv")
    summary = next(r for r in concordance if r["record_type"] == "summary")
    if "36_shared_called_SNPs" not in summary["interpretation"]:
        raise AssertionError("mitochondrial concordance summary changed")


def main() -> None:
    validate_tkgwv2()
    validate_readv2()
    validate_king_ibs0()
    validate_mtdna()
    print("PASS  TKGWV2 16-state table")
    print("PASS  READv2 210-pair cohort table")
    print("PASS  KING/IBS0 genome-wide and chromosome-1 summaries")
    print("PASS  mitochondrial classification and concordance tables")


if __name__ == "__main__":
    main()
