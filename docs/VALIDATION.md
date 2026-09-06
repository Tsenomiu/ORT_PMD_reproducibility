# Validation

Run from the repository root:

```bash
ORT_PCA_VALIDATION_MODE=ci-portable bash reproduce_all.sh
```

Outputs go to `_build/`. For syntax, excluded-file, private-path, credential and
checksum checks alone, run `bash check_repository.sh`. `MANIFEST.sha256` covers
every release file except itself; runtime outputs are excluded.

## Included checks

| Component | Verification |
|---|---|
| DeltaALT | 26 LOCO state intervals, 90 primary paired contrasts and 24 common-callable sensitivity contrasts; synthetic block-jackknife tests |
| Concordance | Dosage/best-guess column identity, dosage reconstruction and byte-exact chromosome-jackknife outputs |
| Mitochondrial | 12 tests, invalid-input rejection and corrected 36/37/36 SNP counts; target-site and rCRS-spacer summaries |
| TKGWV2 | 16 reported-versus-reproduced comparisons, including HRC 0.2356 and 4,517,214 SNPs |
| READv2 | 43 checks of primary and sensitivity summaries and normalization arithmetic |
| Table 1 | 16 generated rows and 80 cell-to-source mappings; comparison with a separately supplied Word manuscript is optional |
| PCA | Numerical tolerance 1e-12, structural checks and three raster comparisons |
| Main outputs | Eight source/status checks covering Figures 1–7 and Table 1; Figure 1 is external PowerPoint artwork |

The default `canonical-strict` PCA mode requires byte-identical raster checksums
from the archived Poppler 25.06.0 / TeX Live 2024 environment. `ci-portable`
retains all 17 numerical/structural checks, verifies the canonical fixtures and
checks raster format, dimensions and visual-equivalence thresholds. Both modes
require `pdftex` and `pdftoppm` for the raster checks. See the
[PCA validation instructions](../workflows/10_pca/recovery/README.md).

## Provenance and limits

These public checks verify compact derived data, not a complete raw-read rerun.
The [workflow map](workflow_manifest.md) distinguishes recovered and reconstructed
code. [Reproducibility limitations](reproducibility_limitations.md) describe
external inputs and the unpreserved pseudo-haploid random states.

The corrected mitochondrial external run retained 45,268 ORT15 and 40,434 ORT16
MT alignments, with no exact repeated SAM records. Corrected VCFs reproduced the
36/37 caller-emitted SNP counts and 36 exact shared calls; both assignments remain
D4o1. The [output checksums](../workflows/09_mtdna/corrected_output_checksums.tsv)
identify those external files. C7028T remains unconfirmed single-read evidence;
no coverage at A4215G is not reference evidence, and `MT:3106 CN>C` is excluded
as an artificial rCRS-spacer record.

Release changes are recorded in the [changelog](../CHANGELOG.md). Historical
reports and per-release inventories remain unchanged in tags
`v1.1.0`, `v1.1.1` and `v1.2.0`.
