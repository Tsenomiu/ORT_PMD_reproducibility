# Compact derived inputs

These distributable summaries support the article figures and numerical checks.
They are derived results, not raw sequence reads or individual-level reference-panel
genotypes.

- `damage/alt_fraction_recomputed_dp3_adsum.csv`: the 26 correction-state rows used
  by main Figure 2; `alt_fraction_summary.tsv` is the equivalent compact table used
  by Supplementary Table S2.
- `imputation/`: the 16-row treatment summary, per-chromosome concordance counts and
  expected paired-jackknife outputs used by main Figure 3 and Table 1.
- `ancibd/`: asymmetric and 64-state IBD1 summaries used by main Figure 4 and
  Supplementary Figures S5--S6.
- `kinship/tkgwv2_pair_results.tsv`: all 16 TKGWV2 correction-state comparisons in
  Supplementary Table S6.
- `kinship/readv2_cemetery_210_pairs.tsv`: the complete 210-pair READv2 output used
  for the cohort counts, ranges, and ORT15–ORT16 parent–offspring result.
- `kinship/readv2_ort15_ort16.tsv`: the primary and autosome-only ORT15–ORT16 rows.
- `kinship/ibs0_positive_control.tsv`: the same-individual chromosome-1 nested
  read-set control.
- `mtdna/`: HaploGrep classifications and the compact ORT15--ORT16 call comparison.
- `pca/`: query coordinates, eigenvalues, query manifest, population-level reference
  aggregates and expected PCA metrics.
- `qc/`: library QC, read-set comparison, radiocarbon and contamination summaries.

Large genomic intermediates, raw images and third-party reference resources are
intentionally absent.
