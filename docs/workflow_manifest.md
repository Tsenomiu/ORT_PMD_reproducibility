# Recovered workflow manifest

This document identifies the publication copies curated for the proposed v1.1.0
release. Private provenance originals remain outside the repository. Machine-
specific paths, server addresses, licensed individual-level reference
coordinates, and large genomic inputs are not included.

| Workflow | Public entry point | Provenance classes | Validation boundary |
|---|---|---|---|
| TKGWV2 | `workflows/08_kinship/tkgwv2/` | Original recovered code, exact upstream compatibility patch, runner reconstructed from successful logs, processed source summaries, generated validation output, documentation | All 16 reported rows reproduced exactly from retained processed PED/MAP inputs; BAM-to-PED randomness is outside the deterministic boundary |
| READv2 | `workflows/08_kinship/readv2/` and `data/processed/readv2/` | Original recovered downstream commands, reconstructed/path-neutral commands, processed normalization data, generated validation output, documentation | Reported primary, autosome-only, and transversion-only arithmetic validated; exact primary BAM-to-genotype command/random state not recovered |
| PCA | `workflows/10_pca/recovery/` and `data/processed/pca/` | Path-neutral preparation/projection records, original parameter settings, processed query-only coordinates and aggregates, generated validation output, documentation | 20 checks passed; Figure 7 and Supplementary Figures S7–S8 raster comparisons matched exactly |
| Table 1 | `workflows/06_concordance/table1/`, `data/processed/table1/`, and `docs/validation/` | Original recovered generator, existing source summaries, generated output, cell-level validation, documentation | 16/16 generated rows and 80/80 Word table cells validated |

The source-data and script mapping for all seven main figures and Table 1 is in
[`main_output_source_manifest.tsv`](main_output_source_manifest.tsv).

## Classification definitions

- **Original recovered code:** an unchanged publication copy of an executed or
  historically retained author file.
- **Reconstructed from logs:** a path-neutral runner assembled from preserved
  successful commands; it is not presented as the original wrapper.
- **Generated validation output:** a newly generated comparison or report used
  to test the curated package.
- **Processed source data:** compact derived values or manifests that support
  verification without redistributing raw or licensed genomic data.
- **Documentation:** explanatory provenance, instructions, limitations, or
  release-preparation records.
