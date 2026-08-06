# Figure 5 — imputation concordance

The generator uses per-MAF-bin dosage summaries, the treatment summary table, and an
independent GCsV harvest. It excludes MAF bin 0 and checks all 16 plotted NRD values
before writing the PDF.

```bash
python figures/figure_05_imputation/make_figure.py \
  --summary data/summary/imputation/imputation_summary.csv \
  --output _build/figures/figure_05_imputation.pdf
```
