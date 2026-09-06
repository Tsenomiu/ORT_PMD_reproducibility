# Validation records

| File | Classification | Purpose |
|---|---|---|
| `table1_cell_provenance.tsv` | generated validation output | Maps all 80 Table 1 data cells to the included imputation summary |
| `table1_validation_results.tsv` | generated validation output | Records the 16-row and 80-cell checks |
| `validate_main_output_sources.py` | generated validation output | Checks repository source files and scripts for Figures 2–7 and Table 1, plus the documented external author-artwork status of Figure 1 |

The protected Word manuscript is not stored in this repository. The optional
Table 1 Word-cell check is documented in
`workflows/06_concordance/table1/README.md` and reads the manuscript without
modifying it.
