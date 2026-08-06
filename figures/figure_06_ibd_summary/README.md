# Figure 6 and Supplementary Figure S4 — IBD summaries

The script validates the asymmetric subset, all 64 cross-treatment rows, the callable
denominator, strict length thresholds, and the 220-SNP/cM density filter.

```bash
python figures/figure_06_ibd_summary/make_figures.py \
  --asymmetric data/summary/ancibd/asymmetric_fulludg_raw.tsv \
  --summary data/summary/ancibd/cross_treatment_pair_summary.tsv \
  --resource-manifest figures/_shared/ancibd_v62/resource_manifest_minimal.json \
  --output-dir _build/figures
```
