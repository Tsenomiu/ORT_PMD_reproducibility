#!/usr/bin/env python3
"""Small regression test guarding GLIMPSE best-guess/dosage column order."""
from __future__ import annotations

import tempfile
from pathlib import Path

from summarize_concordance import weighted_r2


with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "example.rsquare.grp.txt"
    path.write_text("0 10 0.0005 0.10 0.20\n1 2 0.005 0.30 0.50\n2 3 0.015 0.40 0.70\n")
    dosage, best, count = weighted_r2(path)
    assert count == 5
    assert abs(dosage - 0.62) < 1e-12
    assert abs(best - 0.36) < 1e-12

