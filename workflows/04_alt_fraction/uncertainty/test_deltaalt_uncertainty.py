#!/usr/bin/env python3
"""Fast schema, critical-value and Figure 2 input-guard tests."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


primary = load("compute_deltaalt_uncertainty", HERE / "compute_deltaalt_uncertainty.py")
common = load("compute_deltaalt_common_callable", HERE / "compute_deltaalt_common_callable.py")
figure = load("figure02_deltaalt", ROOT / "figures/figure_02_alt_fraction/make_figure.py")


class SourceSchemaTests(unittest.TestCase):
    def source(self, directory: Path, public: bool) -> Path:
        path = directory / ("public.tsv" if public else "round45.tsv")
        if public:
            header = "CHROM\tPOS\tREF\tALT\tDP\tAD_REF\tAD_ALT"
            rows = [f"{chrom}\t1\tC\tT\t3\t2\t1" for chrom in range(1, 23)]
        else:
            header = "CHR\tPOS\tREF\tALT\tDP\tREF_COUNT\tALT_COUNT\tALT_FRAC\tCLASS"
            rows = [
                f"{chrom}\t1\tC\tT\t3\t2\t1\t0.333333\tdamage_transition"
                for chrom in range(1, 23)
            ]
        path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
        return path

    def test_primary_accepts_round45_and_public_headers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for public in (False, True):
                with self.subTest(public=public):
                    result = primary.stream_one(
                        (("ORT15", "fu", "raw"), str(self.source(directory, public)))
                    )
                    self.assertEqual(result["data_rows"], 22)
                    self.assertEqual(result["direct"][1], 22)

    def test_common_reader_accepts_public_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            reader = common.CheckedReader(
                ("ORT15", "fu", "raw"), self.source(Path(tmp), True)
            )
            result = reader.finish()
            self.assertEqual(result["data_rows"], 22)
            self.assertEqual(result["unique_variant_keys"], 22)

    def test_fixed_student_t_constants(self) -> None:
        self.assertEqual(primary.t_critical_975(21), 2.07961384473)
        self.assertEqual(primary.t_critical_975(576), 1.96409103161)
        self.assertEqual(primary.t_critical_975(288), 1.96823517361)
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            primary.t_critical_975(20)


class FigureInputGuardTests(unittest.TestCase):
    def test_duplicate_point_and_interval_keys_are_rejected(self) -> None:
        point_source = ROOT / "data/summary/damage/alt_fraction_recomputed_dp3_adsum.csv"
        interval_source = ROOT / "data/summary/damage/deltaalt_uncertainty/state_ci_loco.csv"
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            point = directory / "point.csv"
            interval = directory / "interval.csv"
            point_text = point_source.read_text(encoding="utf-8")
            interval_text = interval_source.read_text(encoding="utf-8")
            point.write_text(point_text + point_text.splitlines()[1] + "\n", encoding="utf-8")
            interval.write_text(
                interval_text + interval_text.splitlines()[1] + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "Duplicate state"):
                figure.read_rows(point)
            with self.assertRaisesRegex(ValueError, "Duplicate state"):
                figure.read_intervals(interval)


if __name__ == "__main__":
    unittest.main()
