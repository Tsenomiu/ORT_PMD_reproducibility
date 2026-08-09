# Supplementary Figure S4 — PMD profiles after correction

Included aggregate mapDamage tables cover uncorrected, Trim-5, Rescale-5, and
bamRefine-5 treatments for both individuals and library types. Rescale-5 reuses the
uncorrected mismatch table because rescaling changes base qualities rather than the
observed mismatch counts.

The full-UDG rows use merged pre-deduplication collapsed-BAM tables. The non-UDG
rows use post-deduplication collapsed-BAM tables, as labelled in the figure. These
profiles support within-row comparisons among correction states; direct quantitative
comparison between full-UDG and non-UDG rows would also reflect the different
duplicate-processing stages.

```bash
python figures/supplementary_04_pmd_corrected/make_figure.py \
  --output _build/figures/supplementary_04_pmd_corrected.pdf
```
