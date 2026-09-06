#!/usr/bin/env python3
"""Unit tests for the corrected public mitochondrial workflow."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load_module(name: str, filename: str):
    specification = importlib.util.spec_from_file_location(name, HERE / filename)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {filename}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


manifest_module = load_module("mtdna_validate_manifest", "validate_manifest.py")
support_module = load_module("mtdna_target_support", "summarize_target_sites.py")


class ManifestValidationTests(unittest.TestCase):
    def write_manifest(self, directory: Path, rows: list[list[str]], header: str | None = None) -> Path:
        path = directory / "manifest.tsv"
        lines = [header or "sample\tread_set_state\tinput_scope\tbam"]
        lines.extend("\t".join(row) for row in rows)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def valid_rows(self) -> list[list[str]]:
        return [
            ["ORT15", "full_udg_collapsed_only", "combined_libraries_a_b_c", "/input/ORT15_fu.1.bam"],
            ["ORT16", "full_udg_collapsed_only", "combined_libraries_a_b_c", "/input/ORT16_fu.1.bam"],
        ]

    def test_valid_canonical_state_one(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = manifest_module.validate_manifest(self.write_manifest(Path(tmp), self.valid_rows()))
        self.assertEqual([row["sample"] for row in rows], ["ORT15", "ORT16"])

    def test_rejects_nested_state_filenames(self) -> None:
        for nested_name in (
            "ORT15_fu.2.bam", "ORT15_fu.3.bam",
            "ORT15_fu.2.sorted.bam", "ORT15_fu.3-reheader.bam",
        ):
            with self.subTest(name=nested_name), tempfile.TemporaryDirectory() as tmp:
                rows = self.valid_rows()
                rows[0][-1] = f"/input/{nested_name}"
                with self.assertRaisesRegex(ValueError, "nested-state"):
                    manifest_module.validate_manifest(self.write_manifest(Path(tmp), rows))

    def test_rejects_renamed_bam_with_wrong_checksum(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            renamed = directory / "ORT15.canonical-looking.bam"
            renamed.write_bytes(b"not the verified collapsed-only BAM")
            rows = self.valid_rows()[:1]
            rows[0][-1] = str(renamed)
            with self.assertRaisesRegex(ValueError, "BAM SHA-256"):
                manifest_module.validate_manifest(
                    self.write_manifest(directory, rows), require_files=True
                )

    def test_rejects_noncollapsed_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.valid_rows()
            rows[0][1] = "collapsed_plus_truncated"
            with self.assertRaisesRegex(ValueError, "only 'full_udg_collapsed_only'"):
                manifest_module.validate_manifest(self.write_manifest(Path(tmp), rows))

    def test_rejects_duplicate_sample_or_bam(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.valid_rows()
            rows[1][0] = "ORT15"
            with self.assertRaisesRegex(ValueError, "duplicate sample"):
                manifest_module.validate_manifest(self.write_manifest(Path(tmp), rows))
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.valid_rows()
            rows[1][-1] = rows[0][-1]
            with self.assertRaisesRegex(ValueError, "BAM reused"):
                manifest_module.validate_manifest(self.write_manifest(Path(tmp), rows))

    def test_rejects_wrong_header_and_missing_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_manifest(Path(tmp), self.valid_rows(), header="sample\tbam")
            with self.assertRaisesRegex(ValueError, "expected manifest header"):
                manifest_module.validate_manifest(path)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "does not exist"):
                manifest_module.validate_manifest(
                    self.write_manifest(Path(tmp), self.valid_rows()), require_files=True
                )

    def test_rejects_extra_data_field(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_manifest(Path(tmp), self.valid_rows())
            text = path.read_text(encoding="utf-8").replace(
                "/input/ORT15_fu.1.bam\n", "/input/ORT15_fu.1.bam\textra\n"
            )
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "too many"):
                manifest_module.validate_manifest(path)


class TargetSupportTests(unittest.TestCase):
    def test_cigar_projection(self) -> None:
        self.assertEqual(
            support_module.base_at_ref_pos("ACGT", "IIII", 100, "4M", 102),
            ("G", 40),
        )
        self.assertEqual(
            support_module.base_at_ref_pos("AACGT", "IIIII", 100, "2M1I2M", 103),
            ("T", 40),
        )
        self.assertEqual(
            support_module.base_at_ref_pos("ACGT", "IIII", 100, "2M1D2M", 102),
            ("*", None),
        )

    def test_public_output_has_no_read_name_field(self) -> None:
        self.assertNotIn("qname", support_module.OUTPUT_FIELDS)
        self.assertNotIn("bam", support_module.OUTPUT_FIELDS)

    def test_sam_flag_quality_strand_and_distinct_name_counting(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            fake = directory / "samtools"
            sam_rows = [
                "molA\t0\tMT\t7028\t37\t1M\t*\t0\t0\tT\tI",
                "molA\t16\tMT\t7028\t37\t1M\t*\t0\t0\tT\tI",
                "ref1\t0\tMT\t7028\t35\t1M\t*\t0\t0\tC\tH",
                "other1\t0\tMT\t7028\t40\t1M\t*\t0\t0\tA\tI",
                "dup\t1024\tMT\t7028\t60\t1M\t*\t0\t0\tT\tI",
                "secondary\t256\tMT\t7028\t60\t1M\t*\t0\t0\tT\tI",
                "low_bq\t0\tMT\t7028\t60\t1M\t*\t0\t0\tT\t!",
            ]
            fake.write_text(
                "#!/bin/sh\nprintf '%b\\n' "
                + " ".join(repr(row) for row in sam_rows)
                + "\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            row = support_module.summarize_site(
                str(fake), directory / "input.bam", "ORT16", "MT", 7028,
                "C", "T", 30, 30,
            )
            self.assertEqual(
                tuple(row[field] for field in (
                    "qualifying_ref_reads", "qualifying_alt_reads",
                    "qualifying_other_reads", "distinct_alt_read_names",
                    "alt_forward", "alt_reverse",
                )),
                (1, 2, 1, 1, 1, 1),
            )
            self.assertAlmostEqual(float(row["alt_fraction_ref_alt"]), 2 / 3)
            self.assertEqual((row["alt_mapq_min"], row["alt_mapq_max"]), (37, 37))
            self.assertEqual((row["alt_baseq_min"], row["alt_baseq_max"]), (40, 40))


class CallSummaryTests(unittest.TestCase):
    @staticmethod
    def write_vcf(path: Path, sample: str, include_7028: bool) -> None:
        rows = []
        for index in range(36):
            position = 1000 + index
            rows.append(f"MT\t{position}\t.\tA\tG\t50\t.\tDP=10\tGT\t1")
        if include_7028:
            rows.append("MT\t7028\t.\tC\tT\t8.13869\t.\tDP=1;DP4=0,0,1,0\tGT\t1")
        if sample == "ORT15":
            spacer = "DP=152;IDV=126;IMF=0.828947;DP4=7,10,65,57"
        else:
            spacer = "DP=148;IDV=131;IMF=0.885135;DP4=6,4,65,63"
        rows.append(f"MT\t3106\t.\tCN\tC\t200\t.\t{spacer}\tGT\t1")
        text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
        path.write_text(text + "\n".join(rows) + "\n", encoding="utf-8")

    def test_corrected_summary_rebuilds_public_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            ort15 = directory / "ORT15.vcf"
            ort16 = directory / "ORT16.vcf"
            output = directory / "output"
            self.write_vcf(ort15, "ORT15", include_7028=False)
            self.write_vcf(ort16, "ORT16", include_7028=True)
            subprocess.run(
                [
                    sys.executable,
                    str(HERE / "summarize_calls.py"),
                    "--ort15-vcf",
                    str(ort15),
                    "--ort16-vcf",
                    str(ort16),
                    "--target-support",
                    str(ROOT / "data/summary/mtdna/mtdna_target_site_support.tsv"),
                    "--output-dir",
                    str(output),
                ],
                check=True,
            )
            self.assertEqual(
                (output / "mtdna_pair_concordance.tsv").read_bytes(),
                (ROOT / "data/summary/mtdna/mtdna_pair_concordance.tsv").read_bytes(),
            )
            self.assertEqual(
                (output / "mtdna_3106_spacer_audit.tsv").read_bytes(),
                (ROOT / "data/summary/mtdna/mtdna_3106_spacer_audit.tsv").read_bytes(),
            )
            with (output / "ORT15_ORT16.shared_exact_snps.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 36)

    def test_rejects_incorrect_7028_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            ort15 = directory / "ORT15.vcf"
            ort16 = directory / "ORT16.vcf"
            support = directory / "support.tsv"
            self.write_vcf(ort15, "ORT15", include_7028=False)
            self.write_vcf(ort16, "ORT16", include_7028=True)
            text = (ROOT / "data/summary/mtdna/mtdna_target_site_support.tsv").read_text(
                encoding="utf-8"
            )
            support.write_text(text.replace("\t1\t0\t1.0\t37", "\t1\t1\t1.0\t37"), encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(HERE / "summarize_calls.py"),
                    "--ort15-vcf",
                    str(ort15),
                    "--ort16-vcf",
                    str(ort16),
                    "--target-support",
                    str(support),
                    "--output-dir",
                    str(directory / "output"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("unexpected aggregate support", completed.stderr)


if __name__ == "__main__":
    unittest.main()
