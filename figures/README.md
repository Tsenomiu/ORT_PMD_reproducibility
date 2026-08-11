# Figure generation

Run all figure generators supported by the included compact inputs from the
repository root:

```bash
python -m pip install -r requirements.txt
bash figures/reproduce_figures.sh _build/figures
```

The script reports explicit skips for Supplementary Figure S2 (controlled
archaeological spatial inputs) and Supplementary Figure S3 (controlled CT volumes
and selected slice coordinates). PCA plots use non-identifying population-level
reference aggregates; exact individual-level reference backgrounds require licensed
AADR data.

## Figure map

| Figure | Subject | Generator | Input status |
|---|---|---|---|
| 1 | Study workflow | `figure_01_workflow/` | Included |
| 2 | Damage-associated allele support and site retention | `figure_02_alt_fraction/` | Included |
| 3 | Cross-library imputation agreement | `figure_03_imputation/` | Included |
| 4 | IBD1 length summaries | `figure_04_ibd_summary/` | Included |
| 5 | IBD1 karyograms | `figure_05_ibd_karyogram/` | Included |
| 6 | Kinship summary | `figure_06_kinship/` | Included |
| 7 | PCA across correction states | `figure_07_pca/` | Aggregate reference background |
| S1 | PMD profiles before correction | `supplementary_01_pmd_uncorrected/` | Included |
| S2 | Geographic location and site plan | Not distributed | Controlled input required |
| S3 | Petrous micro-CT documentation | `supplementary_03_ct/` | Controlled input required |
| S4 | PMD profiles after correction | `supplementary_04_pmd_corrected/` | Included |
| S5 | Cross-treatment IBD1 matrix | `supplementary_05_ibd_matrix/` | Included |
| S6 | IBD1 length-fragmentation diagnostic | `figure_04_ibd_summary/` | Included |
| S7–S8 | Coverage-titration PCA | `supplementary_07_08_pca_titration/` | Aggregate reference background |

## Interpretation and provenance notes

- **Figure 1:** The editable SVG is rendered through the bundled Matplotlib code so
  its embedded typeface matches the data figures.
- **Figure 2:** Alternate-read fractions use combined REF+ALT depth of at least 3;
  retention uses depth of at least 1. Covered sites are treatment-specific, so the
  contrast is descriptive and is not a fixed-common-site test of variant loss.
- **Figure 3:** The uncorrected full-UDG imputed dataset is the fixed comparator for
  every non-UDG correction state from the same individual. It is not known genotype
  truth. The generator checks the 32 weighted concordance summaries, 16 NRD values,
  and 128 plotted MAF-bin rows.
- **Figures 4–5 and S5–S6:** IBD1 displays use strict length and density filters
  documented in the scripts and included resource manifest. The Figure 5 karyogram
  uses segments longer than 12 cM, density above 220 SNP/cM, and no external gap
  merging.
- **Figure 6:** TKGWV2, KING/IBS0, and READv2 use complementary input constructions.
  The READv2 point comes from the primary cemetery-wide 1240k callset; the autosome-
  only sensitivity run returned the same parent–offspring subtype.
- **Figure 7:** The included query coordinates cover 24 imputed and 24 pseudo-haploid
  projections. `inputs/reference_aggregates.npz` contains only non-identifying grid
  summaries and cannot reconstruct individual AADR reference genotypes.
- **Figure S1:** Four aggregate mapDamage tables provide 5-prime C-to-T and 3-prime
  G-to-A frequencies over terminal positions 1–25.
- **Figure S3:** `SCAN_METADATA.md` records acquisition and display settings. The
  exporter selects recorded slices, applies the shared display window, averages five
  adjacent slices, and adds scale bars; it does not alter the depicted anatomy.
- **Figure S4:** Full-UDG rows use merged pre-deduplication collapsed-BAM profiles;
  non-UDG rows use post-deduplication profiles. Quantitative comparisons should be
  made within rows. Rescale-5 reuses the uncorrected mismatch profile because
  rescaling changes base qualities, not observed mismatches.
- **Figures S7–S8:** The included file contains 1,694 ORT titration projections and
  non-identifying aggregate reference backgrounds. `pdftex` is required for the
  deterministic PDF wrapper.

Each generator supports `--help` or documents its required output environment in the
source header. `reproduce_figures.sh` is the recommended entry point because it also
validates the expected row counts, filters, and output names.
