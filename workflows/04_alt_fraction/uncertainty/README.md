# DeltaALT uncertainty (Round45)

This directory preserves the scripts used for the Reduce18F DeltaALT uncertainty analysis. DeltaALT is the mean site-level ALT fraction at C>T/G>A SNPs minus the corresponding mean at transversion SNPs, in percentage points, among sites with `REF+ALT >= 3` reads.

## Analysis design

- Primary estimates use each correction state's own eligible-site set.
- Uncertainty is calculated with paired leave-one-autosome-out (LOCO), 5-Mb, and 10-Mb genomic-block jackknives.
- Paired contrasts omit the same physical block from every term.
- The common-callable sensitivity analysis intersects exact `CHR/POS/REF/ALT` keys and requires the depth threshold in every state contributing to a comparison.
- Exact identical adjacent source rows are preserved for the primary canonical calculation and reported separately; the common-callable intersection collapses exact duplicate keys and rejects conflicting duplicates.

## Included and external inputs

The compact release tables are in `data/summary/damage/deltaalt_uncertainty/`. They include all 26 LOCO state intervals, all 90 primary paired-contrast summaries, all 24 common-callable sensitivity summaries, state summaries, diagnostics, canonical validation, and path-neutral source checksums. No per-site read-support TSV or private server path is included.

The two compute scripts require the external per-site TSV root because those files contain tens of millions of rows. The expected layout is encoded by `source_path()` in `compute_deltaalt_uncertainty.py`; verify external files against `source_checksums.tsv` before a full rerun.

```sh
python3 workflows/04_alt_fraction/uncertainty/compute_deltaalt_uncertainty.py \
  NEW_PRIMARY_OUTPUT --source-root EXTERNAL_GLIMPSE_ROOT --workers 8

python3 workflows/04_alt_fraction/uncertainty/compute_deltaalt_common_callable.py \
  NEW_COMMON_OUTPUT --source-root EXTERNAL_GLIMPSE_ROOT --workers 8

python3 workflows/04_alt_fraction/uncertainty/summarize_for_manuscript.py \
  NEW_PRIMARY_OUTPUT NEW_MANUSCRIPT_SUMMARY
```

The manuscript summarizer requires the complete primary output directory, including delete-one replicates; those large intermediate tables are intentionally outside the compact release. This is an explicit external-input boundary, not an assertion that the full scan can be reproduced from summary CSVs alone.

## Fast public validation

```sh
python3 workflows/04_alt_fraction/uncertainty/validate_deltaalt_uncertainty.py
```

This checks state/contrast completeness, canonical point-estimate agreement, exact headline intervals, the three LOCO state intervals that cross zero, correction-versus-uncorrected interval direction under all three block schemes, and the common-callable sensitivity conclusions.
