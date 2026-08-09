#!/usr/bin/env python3
"""Regression test for posterior-dosage reconstruction and executed precision."""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent

with tempfile.TemporaryDirectory() as directory:
    vcf = Path(directory) / "example.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t1\t.\tA\tG\t.\tPASS\t.\tGP\t0.123456,0.234567,0.641977\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [sys.executable, str(HERE / "add_ds_from_gp.py"), str(vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    data_line = next(line for line in completed.stdout.splitlines()
                     if not line.startswith("#"))
    assert data_line.endswith("0.123456,0.234567,0.641977:1.5185")
    assert '##FORMAT=<ID=DS,Number=1,Type=Float' in completed.stdout
