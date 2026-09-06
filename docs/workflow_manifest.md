# Workflow provenance

The public entry points below distinguish recovered code from reconstructed or
corrected implementations. Private originals, logs and large inputs remain
outside the repository.

| Workflow | Entry point and provenance |
|---|---|
| TKGWV2 | [Package](../workflows/08_kinship/tkgwv2/): recovered compatibility patch and dyads; runner reconstructed from successful logs |
| READv2 | [Package](../workflows/08_kinship/readv2/): recovered downstream commands and path-neutral wrappers; primary BAM-to-genotype command/random state not preserved |
| DeltaALT | [Uncertainty](../workflows/04_alt_fraction/uncertainty/): path-neutral primary/common-callable scripts and source checksums; full per-site tables external |
| Mitochondrial | [Calling](../workflows/09_mtdna/): corrected collapsed-only runner, nested-state guard, exact-call and aggregate-support code; BAM/VCF inputs external |
| PCA | [Recovery](../workflows/10_pca/recovery/): path-neutral preparation/projection records and original settings; public checks use query coordinates and aggregate references |
| Table 1 | [Generator](../workflows/06_concordance/table1/): recovered original generator; separately created cell-level validator |

For check coverage, see [validation](VALIDATION.md). The
[main-output source map](main_output_source_manifest.tsv) identifies figure/table
inputs and generators. Scientific limits are in
[reproducibility limitations](reproducibility_limitations.md).
