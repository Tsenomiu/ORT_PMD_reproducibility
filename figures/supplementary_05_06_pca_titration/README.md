# Supplementary Figures S5 and S6 — coverage-titration PCA

The generator reads 1,694 ORT titration-query coordinates and population-level
aggregate reference backgrounds. Exact individual-level AADR backgrounds require
separately licensed inputs. `pdftex` is required for the deterministic PDF wrapper.

```bash
ORT_FIGURE_OUTPUT_DIR=_build/figures/pca_titration \
  python figures/supplementary_05_06_pca_titration/make_figures.py
```
