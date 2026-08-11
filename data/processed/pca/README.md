# PCA processed-data record

This directory contains publication-safe provenance and validation records for
the recovered PCA workflow. The compact numerical and plotting inputs already
present elsewhere in the repository are intentionally not duplicated:

- `data/summary/pca/pca_ort_queries_matched48.evec`: 48 query rows only;
- `data/summary/pca/pca_aadr_matched48.eval`: final eigenvalues;
- `data/summary/pca/query_manifest_matched48.tsv`: 48-query treatment map;
- `data/summary/pca/reference_population_aggregates.tsv`: population aggregates;
- `data/summary/pca/pca_metrics_matched48.csv`: archived method-distance metrics;
- `figures/figure_07_pca/inputs/`: query coordinates and aggregate background;
- `figures/supplementary_07_08_pca_titration/inputs/`: titration query
  coordinates and aggregate background.

No individual-level AADR coordinates or genotypes are included. The recovered
projection produced 985 reference rows and 48 projected query rows, but the public
coordinate file retains only the 48 query rows. Reference backgrounds used for
plotting are disclosure-controlled population aggregates.

`validation_results.tsv` contains 20 passed checks. The first 17 are numerical or
structural checks; the final three record byte-identical 150-dpi raster
comparisons against the manuscript figures. `input_output_manifest.tsv` maps each
input, script and expected output without exposing machine-specific locations.

These TSV files are the archived `canonical-strict` record. They and
`raster_reference_checksums.tsv` remain unchanged. Canonical rasters were
produced with Poppler 25.06.0 and TeX Live 2024. The failing v1.1.0 GitHub
Actions runs used Poppler 24.02.0, which changed raster anti-aliasing and PNG
bytes without changing numerical outputs or scientific content. The separate
runtime `ci-portable` mode still requires all numerical and structural checks,
RGB PNG format, exact dimensions, and documented visual-equivalence thresholds.
