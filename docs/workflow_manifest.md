# Recovered workflow manifest

This document identifies the publication copies introduced in v1.1.0, retained
in v1.1.1, and extended by corrected mitochondrial and DeltaALT uncertainty
workflows in v1.2.0.
Private provenance originals remain outside the repository.
Machine-specific paths, server addresses, licensed individual-level reference
coordinates, and large genomic inputs are not included.

| Workflow | Public entry point | Provenance classes | Validation boundary |
|---|---|---|---|
| TKGWV2 | `workflows/08_kinship/tkgwv2/` | Original recovered code, exact upstream compatibility patch, runner reconstructed from successful logs, processed source summaries, generated validation output, documentation | All 16 reported rows reproduced exactly from retained processed PED/MAP inputs; BAM-to-PED randomness is outside the deterministic boundary |
| READv2 | `workflows/08_kinship/readv2/` and `data/processed/readv2/` | Original recovered downstream commands, reconstructed/path-neutral commands, processed normalization data, generated validation output, documentation | Reported primary, autosome-only, and transversion-only arithmetic validated; exact primary BAM-to-genotype command/random state not recovered |
| DeltaALT uncertainty | `workflows/04_alt_fraction/uncertainty/` and `data/summary/damage/deltaalt_uncertainty/` | Path-neutral Round45 primary and common-callable scripts, compact jackknife summaries, sanitized source checksums, Figure 2 interval input and validator | All 26 LOCO state intervals, 90 primary contrasts and 24 sensitivity contrasts validated; full per-site TSV scan remains an external-data boundary |
| Mitochondrial calling | `workflows/09_mtdna/` and `data/summary/mtdna/` | Corrected path-neutral runner, nested-state guard, aggregate support and exact-call comparison code, processed summaries, checksums and documentation | D4o1 and scores preserved; 36/37 caller-emitted SNPs and 36 exact shared calls validated; BAM/VCF regeneration requires controlled inputs |
| PCA | `workflows/10_pca/recovery/` and `data/processed/pca/` | Path-neutral preparation/projection records, original parameter settings, processed query-only coordinates and aggregates, generated validation output, documentation | All 17 numerical/structural checks remain strict; canonical validation requires three exact raster SHA matches, while CI requires exact format/dimensions and documented visual equivalence |
| Table 1 | `workflows/06_concordance/table1/`, `data/processed/table1/`, and `docs/validation/` | Original recovered generator, existing source summaries, generated output, cell-level validation, documentation | 16/16 generated rows and 80/80 Word table cells validated |

The source-data and script mapping for all seven main figures and Table 1 is in
[`main_output_source_manifest.tsv`](main_output_source_manifest.tsv).

## Classification definitions

- **Original recovered code:** an unchanged publication copy of an executed or
  historically retained author file.
- **Reconstructed from logs:** a path-neutral runner assembled from preserved
  successful commands; it is not presented as the original wrapper.
- **Corrected path-neutral code:** a publication-safe implementation validated
  against the preserved corrected outputs after a documented defect was removed.
- **Generated validation output:** a newly generated comparison or report used
  to test the curated package.
- **Processed source data:** compact derived values or manifests that support
  verification without redistributing raw or licensed genomic data.
- **Documentation:** explanatory provenance, instructions, limitations, or
  release-preparation records.
