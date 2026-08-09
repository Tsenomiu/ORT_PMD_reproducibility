# Figure 6 — kinship summary

`input_values.tsv` contains the plotted TKGWV2, KING/IBS0, and READv2 summaries.
The unrelated IBS0 expectation is the sitewise mean of `2*p^2*q^2`; the full-sibling
guide is one quarter of that value.

The plotted READv2 value comes from the primary cemetery-wide 1240k callset,
which included autosomal and sex-chromosome targets. The ORT15 and ORT16 inputs
were merged full-UDG a--c BAMs with 6 bp soft-clipped from each read end. The
autosome-only ORT6--ORT27 rerun gave the same first-degree parent--offspring
classification. READv2 is therefore a complementary cohort estimate, not a
matched comparison with the TKGWV2 and KING input constructions.

```bash
python figures/figure_06_kinship/make_figure.py \
  --output _build/figures/figure_06_kinship.pdf
```
