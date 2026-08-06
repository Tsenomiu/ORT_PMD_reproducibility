# 06 Fixed-comparator concordance and chromosome jackknife

For each individual, the unmodified full-UDG imputed VCF is held fixed as the
matched comparator for all non-UDG treatments. It is not an independent genotype
truth.

`add_ds_from_gp.py` reconstructs posterior dosage as
`DS = GP(0/1) + 2*GP(1/1)`. GLIMPSE2_concordance then receives target dosage,
comparator genotype, and reference-panel allele frequencies. The analysis uses:

```text
--gt-val --af-tag AF
--bins 0 0.001 0.01 0.02 0.05 0.10 0.20 0.30 0.40 0.50
```

Dosage and best-guess r-squared are count-weighted over bins 1–8 (MAF at least
0.1%); NRD uses all evaluated sites.

`run_per_chromosome_concordance.sh` applies the same comparison separately to
chromosomes 1–22. `compute_jackknife.py` reconstructs the paired leave-one-
chromosome-out NRD comparison from the included chromosome counts. Run the compact
regression check with `bash reproduce_all.sh` from the repository root.
