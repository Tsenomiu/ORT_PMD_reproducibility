# Figure 2 — damage-associated allele support and site retention

The script reads all 26 treatment rows from
`data/summary/damage/alt_fraction_recomputed_dp3_adsum.csv`. Lines show the treatment-specific
mean alternate-read fraction at C→T/G→A SNPs minus the corresponding mean at
transversion SNPs. Pale bars show covered-site retention relative to the uncorrected
library of the same individual and library type.

Fraction means use combined REF+ALT depth of at least 3. Retention uses combined
REF+ALT depth greater than 0. Callable sites are treatment-specific; this figure is
therefore descriptive and is not a fixed-common-site test of genuine variant loss.

```bash
python figures/figure_02_alt_fraction/make_figure.py \
  data/summary/damage/alt_fraction_recomputed_dp3_adsum.csv \
  _build/figures/figure_02_alt_fraction.pdf
```
