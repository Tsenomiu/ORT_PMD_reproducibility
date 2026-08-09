# Figure 7 — matched-48 PCA

Inputs contain 24 imputed and 24 pseudo-haploid ORT query coordinates, PCA
eigenvalues, and population-level aggregate reference backgrounds. Individual AADR
reference coordinates are not distributed, so exact reconstruction of the
individual-level background requires separately licensed AADR inputs.

```bash
ORT_FIGURE_OUTPUT=_build/figures/figure_07_pca.pdf \
  python figures/figure_07_pca/make_figure.py
```
