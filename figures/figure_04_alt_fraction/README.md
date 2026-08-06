# Figure 4 — alternative-allele fraction

The plotted values are read from `data/summary/damage/alt_fraction_summary.tsv`.
Transition means use sites with analysed depth at least 3; retained-site percentages
use all covered sites.

```bash
python figures/figure_04_alt_fraction/make_figure.py \
  --input data/summary/damage/alt_fraction_summary.tsv \
  --output _build/figures/figure_04_alt_fraction.pdf
```
