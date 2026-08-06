# Supplementary Figure S2 — PMD profiles after correction

Included aggregate mapDamage tables cover uncorrected, Trim-5, Rescale-5, and
bamRefine-5 treatments for both individuals and library types. Rescale-5 reuses the
uncorrected mismatch table because rescaling changes base qualities rather than the
observed mismatch counts.

```bash
python figures/supplementary_02_pmd_corrected/make_figure.py \
  --output _build/figures/supplementary_02_pmd_corrected.pdf
```
