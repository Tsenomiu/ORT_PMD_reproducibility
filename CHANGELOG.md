# Changelog

## v1.1.1

- Add an explicit `ci-portable` PCA raster-validation mode for GitHub Actions.
- Retain `canonical-strict` as the default, including all 20 PCA checks and all
  three exact canonical raster SHA-256 comparisons.
- Document the renderer-dependent anti-aliasing difference between canonical
  Poppler 25.06.0 / TeX Live 2024 output and Poppler 24.02.0 in the failing
  v1.1.0 GitHub Actions runs.
- No scientific data, analysis results, parameters, source values, figure
  content, tables, or conclusions changed.

## v1.1.0

- Add the recovered TKGWV2 v1.0b compatibility patch, a path-neutral
  processed-input runner, input and environment manifests, and exact validation
  of all 16 reported comparisons.
- Add READv2 downstream commands, normalization data, and validators for the
  primary, autosome-only, and transversion-only calculations, with the missing
  primary BAM-to-genotype command and random state disclosed.
- Add the recovered TEST/RC PCA preparation, projection, metric, and plotting
  records in publication-safe form, together with 20 validation checks and
  three exact raster comparisons.
- Add the recovered original Table 1 generator, its generated output, and the
  16-row/80-cell provenance check.
- Add main-output source-data, workflow, validation, file, and exclusion
  manifests for the proposed release.

## v1.0.1

- Simplified the reproducibility package while preserving the reported results
  and supported figure regeneration.
- Added compact KING/IBS0 and READv2 summary evidence and improved portability
  checks.
