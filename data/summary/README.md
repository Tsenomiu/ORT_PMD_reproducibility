# Compact derived inputs

These files support the plotting and numerical checks included in this repository.
They are derived summaries, not raw sequence or genotype data.

- `damage/alt_fraction_summary.tsv`: alternative-allele fractions and retained-site
  counts used by Figure 4.
- `imputation/imputation_summary.csv`: treatment-level imputation metrics used by
  Figure 5.
- `imputation/perchr_gcsv_rows.txt`: per-chromosome GLIMPSE concordance counts used
  to reconstruct the paired jackknife.
- `imputation/perchr_nrd_counts.csv` and `jackknife_deltaNRD.csv`: expected outputs
  used by the jackknife regression check.
- `ancibd/asymmetric_fulludg_raw.tsv` and
  `ancibd/cross_treatment_pair_summary.tsv`: compact focal-pair summaries used by
  Figures 6, S3, and S4.
- `pca/`: focal query coordinates, eigenvalues, query manifest, population-level
  reference summaries, and expected PCA metrics used by the PCA verification step.

Large genomic intermediates and third-party reference data are intentionally absent.
