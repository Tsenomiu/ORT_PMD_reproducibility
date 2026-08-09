# Figure 3 — cross-library imputation agreement

The generator uses per-MAF-bin dosage summaries, the treatment summary table, and an
independent GCsV harvest. Panels a–b plot dosage r-squared for MAF bins 1–8
(MAF at least 0.1%). Panels c–d plot genome-wide NRD from all concordance-evaluable
non-reference genotype comparisons across bins 0–8, with axes beginning at zero.
The same uncorrected full-UDG imputed dataset is the comparator for every non-UDG
treatment from an individual; it is not known genotype truth. The script recomputes
and checks all 32 count-weighted dosage/best-guess summary values, checks all 16 NRD
values, and exports all 128 per-bin plotting rows before writing the PDF.

```bash
python figures/figure_03_imputation/make_figure.py \
  --summary data/summary/imputation/imputation_summary.csv \
  --output _build/figures/figure_03_imputation.pdf
```
