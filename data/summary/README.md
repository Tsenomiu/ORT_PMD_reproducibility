# Compact derived inputs

These distributable summaries support the article figures and numerical checks.
They are derived results, not raw sequence reads or individual-level reference-panel
genotypes.

- `damage/alt_fraction_recomputed_dp3_adsum.csv`: the 26 correction-state rows used
  by main Figure 2 and Supplementary Table S2.
- `damage/deltaalt_uncertainty/`: compact Round45 LOCO state intervals, primary
  paired jackknife contrasts, common-callable sensitivity summaries, diagnostics,
  canonical validation and sanitized source checksums used by Reduce18F/Figure 2.
- `imputation/`: the 16-row treatment summary, per-chromosome concordance counts and
  expected paired-jackknife outputs used by main Figure 3 and Table 1.
- `ancibd/`: asymmetric and 64-state IBD1 summaries used by main Figure 4 and
  Supplementary Figures S5--S6.
- `kinship/tkgwv2_pair_results.tsv`: all 16 TKGWV2 correction-state comparisons in
  Supplementary Table S6.
- `kinship/readv2_cemetery_210_pairs.tsv`: the complete 210-pair READv2 output used
  for the cohort counts, ranges, and ORT15–ORT16 parent–offspring result.
- `kinship/readv2_ort15_ort16.tsv`: the primary and autosome-only ORT15–ORT16 rows.
- `kinship/king_ibs0_summary.tsv`: genome-wide and high-confidence chromosome-1
  KING/IBS0 summaries.
- `kinship/ibs0_positive_control.tsv`: the same-individual chromosome-1 nested
  read-set control.
- `mtdna/`: corrected collapsed-only HaploGrep classifications, the compact
  ORT15--ORT16 exact-call comparison, aggregate 4215/7028 support and the
  non-biological rCRS-spacer audit. No read names are included.
- `pca/`: query coordinates, eigenvalues, query manifest, population-level reference
  aggregates and expected PCA metrics.
- `qc/`: library QC, read-set comparison, radiocarbon and contamination summaries.

Large genomic intermediates, raw images and third-party reference resources are
intentionally absent.
