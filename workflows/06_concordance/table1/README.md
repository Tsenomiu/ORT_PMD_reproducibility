# Table 1 generation and validation

This directory restores the numerical generator used for manuscript Table 1 and
adds a publication validator. No exporter was invented and represented as
historical code.

## Provenance labels

- `generate_table1.py`: **original recovered code**, copied without changing its
  calculations.
- `validate_table1.py`: **generated validation code** created during the v1.1.0
  recovery pass.
- `../../../figures/figure_03_imputation/impute_gcsv.json` and
  `../../../figures/figure_03_imputation/rsquare_all16.json`: **source data**.
- `../../../data/processed/table1/table1_generator_output.csv`: **generated
  validation output** from the recovered generator.
- `../../../data/summary/imputation/imputation_summary.csv`: **processed source
  data** used by the manuscript and figures.
- `../../../docs/validation/table1_cell_provenance.tsv`: **generated validation
  output**, mapping each of the 80 manuscript cells to the processed summary.

The generator excludes MAF bin 0, computes count-weighted means over bins 1–8,
keeps dosage and best-guess r² separate, and reports 16 sample/treatment rows.
`Raw (uncorrected)` in its direct output is the same presentation row called
`Uncorrected` in the Word table and processed summary.

## Regenerate the source table

From the repository root:

```bash
python3 workflows/06_concordance/table1/generate_table1.py \
  figures/figure_03_imputation/impute_gcsv.json \
  figures/figure_03_imputation/rsquare_all16.json \
  _build/tables/table1.csv
```

Validate the 16 generated rows and the stored 80-cell provenance:

```bash
python3 workflows/06_concordance/table1/validate_table1.py
```

When the protected manuscript is available, run the full read-only Word check in
an environment containing `python-docx`:

```bash
python3 workflows/06_concordance/table1/validate_table1.py \
  --manuscript /path/to/ORT_PMD_reduce16_with_figures_tables.docx
```

Recovery validation passed for **16/16 rows** and **80/80 Word table cells**.
