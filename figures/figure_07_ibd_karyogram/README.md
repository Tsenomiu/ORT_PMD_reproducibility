# Figure 7 — IBD1 karyogram

The karyogram shows density-filtered all-SNP IBD1 intervals strictly longer than
12 cM. The included manifest supplies chromosome map bounds and the callable span.

```bash
python figures/figure_07_ibd_karyogram/make_figure.py \
  --base figures/_shared/ancibd_v62 \
  --resource-manifest figures/_shared/ancibd_v62/resource_manifest_minimal.json \
  --output _build/figures/figure_07_ibd_karyogram.pdf
```
