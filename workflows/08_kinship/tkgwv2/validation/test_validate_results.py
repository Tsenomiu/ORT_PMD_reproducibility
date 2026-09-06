#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Guard the separation of frozen, plotted and reproduced TKGWV2 results."""

import copy
from pathlib import Path
import tempfile
import unittest

import validate_results as validator


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.expected = validator.load_expected_results()
        self.figure = validator.read_tsv(validator.FIGURE_SOURCE_PATH)
        self.observed = validator.read_tsv(validator.RECOVERED_VALIDATION_PATH)

    def test_independent_paths_and_baseline(self):
        paths = {p.resolve() for p in (
            validator.EXPECTED_PATH, validator.FIGURE_SOURCE_PATH,
            validator.RECOVERED_VALIDATION_PATH,
        )}
        self.assertEqual(len(paths), 3)
        validator.validate_figure_source(self.expected, self.figure)
        report = validator.build_report(self.expected, self.observed)
        self.assertEqual(sum(r["status"] == "PASS" for r in report), 16)

    def test_changed_figure_value_is_rejected(self):
        changed = copy.deepcopy(self.figure)
        changed[0]["HRC"] = "0.0000"
        with self.assertRaisesRegex(SystemExit, "figure source differs"):
            validator.validate_figure_source(self.expected, changed)

    def test_changed_reproduced_value_fails_comparison(self):
        changed = copy.deepcopy(self.observed)
        row = next(r for r in changed if r["sample1"] == "ORT15_fu_raw")
        row["HRC"] = "0.0000"
        report = validator.build_report(self.expected, changed)
        self.assertEqual(sum(r["status"] == "FAIL" for r in report), 1)

    def test_missing_reproduced_pair_is_rejected(self):
        with self.assertRaisesRegex(SystemExit, "16 unique pairs"):
            validator.build_report(self.expected, self.observed[:-1])

    def test_frozen_reference_cannot_be_silently_changed(self):
        with tempfile.TemporaryDirectory(prefix="ort-tkgwv2-test-") as directory:
            path = Path(directory) / "altered_expected.tsv"
            content = validator.EXPECTED_PATH.read_bytes().replace(b"0.2096", b"0.2097", 1)
            path.write_bytes(content)
            with self.assertRaisesRegex(SystemExit, "checksum changed"):
                validator.load_expected_results(path)


if __name__ == "__main__":
    unittest.main()
