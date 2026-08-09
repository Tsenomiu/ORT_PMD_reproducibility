# Figure 5 — IBD1 karyogram

The karyogram shows density-filtered, transition-inclusive IBD1 intervals strictly longer than
12 cM. The included manifest supplies chromosome map bounds and the callable span;
`KARYOGRAM_SEGMENTS_DENSITY220.tsv` and `KARYOGRAM_PAIR_SUMMARY.tsv` contain the plotted
interval and summary values. The article plot uses a strict $>12$-cM display
threshold, $>220$ SNP/cM and no external gap merging.

```bash
python figures/figure_05_ibd_karyogram/make_figure.py \
  --base figures/_shared/ancibd_v62 \
  --resource-manifest figures/_shared/ancibd_v62/resource_manifest_minimal.json \
  --output-dir _build/figures/karyogram
```
