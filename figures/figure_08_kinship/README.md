# Figure 8 — kinship summary

`input_values.tsv` contains the plotted TKGWV2, KING/IBS0, and READv2 summaries.
The unrelated IBS0 expectation is the sitewise mean of `2*p^2*q^2`; the full-sibling
guide is one quarter of that value.

```bash
python figures/figure_08_kinship/make_figure.py \
  --output _build/figures/figure_08_kinship.pdf
```
