# Supplementary Figure S3 — cross-treatment IBD1 matrix

The generator checks the ordered 8-by-8 design, strict 12-cM threshold, 220-SNP/cM
density filter, and callable denominator.

```bash
python figures/supplementary_03_ibd_matrix/make_figure.py \
  --summary data/summary/ancibd/cross_treatment_pair_summary.tsv \
  --resource-manifest figures/_shared/ancibd_v62/resource_manifest_minimal.json \
  --output _build/figures/supplementary_03_ibd_matrix.pdf
```
