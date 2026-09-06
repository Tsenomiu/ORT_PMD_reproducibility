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
    if len(calls) != 2 or {row["sample"] for row in calls} != {"ORT15", "ORT16"}:
        raise AssertionError("unexpected mitochondrial classification rows")
    by_sample = {row["sample"]: row for row in calls}
    for sample, score in (("ORT15", 0.8911), ("ORT16", 0.9107)):
        row = by_sample[sample]
        expected = {
            "input_read_set": "full-UDG collapsed-only .1",
            "caller": "bcftools 1.16 haploid",
            "calling_filters": "MAPQ>=30;BQ>=30;BAQ_disabled;no_global_DP_filter",
            "haplogrep_version": "3.3.2",
            "tree": "phylotree-fu-rcrs@1.3",
            "haplogroup": "D4o1",
            "interpretation": "software-reported assignment",
        }
        for field, value in expected.items():
            if row[field] != value:
                raise AssertionError(f"unexpected {sample} mitochondrial {field}")
        close(float(row["quality_score"]), score)

    concordance = rows(ROOT / "data/summary/mtdna/mtdna_pair_concordance.tsv")
    if len(concordance) != 4:
        raise AssertionError("expected four mitochondrial concordance rows")
    keyed = {(row["record_type"], row["position"]): row for row in concordance}
    if len(keyed) != 4:
        raise AssertionError("duplicate mitochondrial concordance key")
    expected_rows = {
        ("summary", "NA"): (
            "36_called_ALT_SNPs",
            "37_called_ALT_SNPs",
            "NA",
            "NA",
            "36 exact shared POS_REF_ALT calls; ORT15 call set is a subset of ORT16",
        ),
        ("artifact", "3106"): (
            "CN>C",
            "CN>C",
            "152",
            "148",
            "left-anchored deletion of artificial rCRS N spacer at 3107; exclude from biological SNP and discordance counts",
        ),
        ("site", "4215"): (
            "no_qualifying_coverage",
            "no_qualifying_coverage",
            "0",
            "0",
            "old ORT15 A>G call does not survive; one truncated-read molecule was counted twice after nested-state merge",
        ),
        ("site", "7028"): (
            "no_qualifying_coverage",
            "C>T_1_read_call",
            "0",
            "1",
            "ORT16 call has one forward read (MAPQ37; BQ40; AF=1/1) and is not a confirmed pairwise difference",
        ),
    }
    if set(keyed) != set(expected_rows):
        raise AssertionError("unexpected mitochondrial concordance row keys")
    for key, expected in expected_rows.items():
        row = keyed[key]
        observed = tuple(
            row[field]
            for field in ("ORT15", "ORT16", "ORT15_depth", "ORT16_depth", "interpretation")
        )
        if observed != expected:
            raise AssertionError(f"unexpected mitochondrial concordance row: {key}")

    spacer = {row["sample"]: row for row in rows(ROOT / "data/summary/mtdna/mtdna_3106_spacer_audit.tsv")}
    if set(spacer) != {"ORT15", "ORT16"}:
        raise AssertionError("unexpected MT:3106 spacer rows")
    expected_spacer = {
        "ORT15": ("152", "126", "0.828947", "65", "57"),
        "ORT16": ("148", "131", "0.885135", "65", "63"),
    }
    for sample, expected in expected_spacer.items():
        row = spacer[sample]
        observed = tuple(row[field] for field in ("DP", "IDV", "IMF", "ALT_forward", "ALT_reverse"))
        if (row["position"], row["REF"], row["ALT"], row["GT"]) != ("3106", "CN", "C", "1"):
            raise AssertionError(f"unexpected {sample} MT:3106 allele")
        if observed != expected or "non-biological" not in row["interpretation"]:
            raise AssertionError(f"unexpected {sample} MT:3106 support")

    support = {
        (row["sample"], row["position"]): row
        for row in rows(ROOT / "data/summary/mtdna/mtdna_target_site_support.tsv")
    }
    if set(support) != {
        ("ORT15", "4215"),
        ("ORT15", "7028"),
        ("ORT16", "4215"),
        ("ORT16", "7028"),
    }:
        raise AssertionError("unexpected mitochondrial target-site rows")
    for key, row in support.items():
        if row["mapq_min"] != "30" or row["baseq_min"] != "30":
            raise AssertionError(f"unexpected target-site filters: {key}")
    ort16_7028 = support[("ORT16", "7028")]
    expected_7028 = ("0", "1", "1", "1", "0", "1.0", "37", "37", "40", "40")
    observed_7028 = tuple(
        ort16_7028[field]
        for field in (
            "qualifying_ref_reads",
            "qualifying_alt_reads",
            "distinct_alt_read_names",
            "alt_forward",
            "alt_reverse",
            "alt_fraction_ref_alt",
            "alt_mapq_min",
            "alt_mapq_max",
            "alt_baseq_min",
            "alt_baseq_max",
        )
    )
    if observed_7028 != expected_7028:
        raise AssertionError("unexpected ORT16 C7028T aggregate support")
    for key, row in support.items():
        if key != ("ORT16", "7028") and any(
            int(row[field]) != 0
            for field in (
                "qualifying_ref_reads",
                "qualifying_alt_reads",
                "qualifying_other_reads",
                "distinct_alt_read_names",
                "alt_forward",
                "alt_reverse",
            )
        ):
            raise AssertionError(f"unexpected qualifying coverage: {key}")
        if key != ("ORT16", "7028") and any(
            row[field] != "NA"
            for field in (
                "alt_fraction_ref_alt", "alt_mapq_min", "alt_mapq_max",
                "alt_baseq_min", "alt_baseq_max",
            )
        ):
            raise AssertionError(f"unexpected zero-coverage ALT metrics: {key}")

    manifest = rows(ROOT / "workflows/09_mtdna/full_udg_bam_manifest.example.tsv")
    if len(manifest) != 2 or {row["sample"] for row in manifest} != {"ORT15", "ORT16"}:
        raise AssertionError("unexpected mitochondrial example manifest rows")
    if any(
        row["read_set_state"] != "full_udg_collapsed_only"
        or row["input_scope"] != "combined_libraries_a_b_c"
        or row["bam"].endswith((".2.bam", ".3.bam"))
        for row in manifest
    ):
        raise AssertionError("nested or noncanonical mitochondrial example input")

    input_checksums = {
        row["sample"]: row
        for row in rows(ROOT / "workflows/09_mtdna/corrected_input_checksums.tsv")
    }
    expected_inputs = {
        "ORT15": "674f891b87a43fe4de33e0ed189ba11eada8c935904593f3ab563eb422a64d65",
        "ORT16": "c614879709587d95fdb2e4b44b34ffef4be89be7f3307a0f93d305d1570970ba",
    }
    if set(input_checksums) != set(expected_inputs):
        raise AssertionError("unexpected corrected mitochondrial input-checksum rows")
    for sample, expected_sha in expected_inputs.items():
        row = input_checksums[sample]
        if (
            row["read_set_state"] != "full_udg_collapsed_only"
            or row["input_scope"] != "combined_libraries_a_b_c"
            or row["sha256"] != expected_sha
        ):
            raise AssertionError(f"unexpected corrected input provenance for {sample}")

    output_checksums = {
        (row["sample"], row["output_type"]): row["sha256"]
        for row in rows(ROOT / "workflows/09_mtdna/corrected_output_checksums.tsv")
    }
    expected_outputs = {
        ("ORT15", "MT_BAM"): "01b43a77adae34d3e6d4ee05033a6b9fe06bbe8abb976f4e189d750087b70f16",
        ("ORT15", "mtDNA_VCF_GZ"): "9e2a9f5d8c7e695dabfa8f7eceb99ecf14a2d8cc41e0cd649c82c247e7a27f9c",
        ("ORT15", "HaploGrep_TSV"): "f0c833230e8ac72fc289e9289cc95e8f2b06ee5be8a391587f29061523bc06c2",
        ("ORT16", "MT_BAM"): "cf1af68204a4cb57a426ceccdccd8a4787825e637a318081de584d9db1ecfae2",
        ("ORT16", "mtDNA_VCF_GZ"): "a1aa1a28e8d0831d2da5d57a9aa4d4809f3775fcff68b646a46d166da5a7c50a",
        ("ORT16", "HaploGrep_TSV"): "ccd8d89ddf73e57fff6b9820f8e0bb4edddf0f580f879629fa2bdba335cbb9f3",
    }
    if output_checksums != expected_outputs:
        raise AssertionError("corrected mitochondrial output checksums changed")

    environment = {
        row["component"]: row
        for row in rows(ROOT / "workflows/09_mtdna/environment.tsv")
    }
    expected_environment = {
        "samtools_htslib": ("1.13", "NA"),
        "bcftools_htslib": ("1.16", "NA"),
        "hs37d5_fasta": ("GRCh37/hs37d5", "65add55817fc9a5de8221caf36e84bf8670f94bef9cff7889093692b271a765d"),
        "hs37d5_fai": ("GRCh37/hs37d5", "1eab7540d4b62ef0b43b50581d027be37c1ee57da9e9cb75957e750641c8630c"),
        "HaploGrep3_archive": ("3.3.2", "957165143ddd7d6ad9f34cc3af56b784fc236081e031c99bea4a574d6ad70959"),
        "HaploGrep3_jar": ("3.3.2", "b7ac239b55b8253dfe03c299f7257e186c73b9505e0e7f8b7a8b14850b10d3f9"),
        "phylotree-fu-rcrs_tree_yaml": ("1.3", "776dc754f466ab4f516e94742b3a79f33cec1a20b8fa582f7e0543fc0151d8d0"),
        "phylotree-fu-rcrs_tree_xml": ("1.3", "07275e7a32f2d9ddb13252f357559eaa3537f677600e127f5b3a1a2c089d45be"),
    }
    if set(environment) != set(expected_environment):
        raise AssertionError("unexpected corrected mitochondrial environment rows")
    for component, expected in expected_environment.items():
        row = environment[component]
        if (row["version_or_identifier"], row["sha256"]) != expected or not row["role"]:
            raise AssertionError(f"unexpected environment provenance for {component}")


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
