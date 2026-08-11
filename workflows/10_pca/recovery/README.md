# Recovered PCA workflow

Status: **recovered and validated**.

This directory is the publication copy of the recovered PCA workflow. Query
preparation was performed in the TEST analysis environment, and marker
harmonisation, reference merging, EIGENSTRAT conversion and the final projection
were performed in the RC compute environment. Machine-specific paths and host
details have been replaced with command-line inputs.

The private recovery archive preserves the recovered files byte-for-byte. The
files here are curated, path-neutral copies or generated validation records:

| File | Classification | Purpose |
|---|---|---|
| `prepare_imputed_queries.sh` | publication wrapper reconstructed from two recovered scripts | Build the 24-query imputed 1240k PLINK dataset. |
| `prepare_pseudohaploid_queries.sh` | path-neutral publication copy of recovered code | Build the 24-query pseudo-haploid 1240k PLINK dataset. |
| `run_recovered_matched48_projection.sh` | path-neutral publication wrapper for the recovered RC run | Preserve the recovered query-BIM sanitisation step, run the shared projection implementation and write query/build records. |
| `query_input_manifest.example.tsv` | documentation template | Records the exact 24 treatment states and user-supplied input paths. |
| `keep_pops_ref.txt` | original recovered parameter data | The 154 HGDP/SGDP populations used to define the reference axes. |
| `smartpca_matched48.par.template` | path-neutral publication copy of the recovered parameter file | Documents the final `smartpca` settings. |
| `environment.tsv` | recovered environment record | Software versions used for query preparation and final projection. |
| `validate_pca_recovery.py` | generated validation code | Recalculates all 20 archived checks. |
| `run_validation_and_figures.sh` | publication wrapper | Recalculates metrics, rebuilds three figures and checks their 150-dpi rasters. |
| `canonical_rasters/` | canonical validation fixtures | Byte-exact manuscript rasters whose SHA-256 values are recorded in `data/processed/pca/raster_reference_checksums.tsv`. |

The projection wrapper calls the existing
`workflows/10_pca/run_matched48_projection.sh`. It uses the existing
`harmonize_markers.py` and `label_eigenstrat.py` helpers, so those files are not
duplicated here. The method-distance calculation is the existing
`workflows/10_pca/verify_pca.py`. The plotting code and compact plotting inputs
remain in `figures/figure_07_pca/` and
`figures/supplementary_07_08_pca_titration/`.

## External inputs

Copy `query_input_manifest.example.tsv` outside the repository and replace the
two path columns with the local 24 imputed VCF and 24 BAM locations. The scripts
also require:

- an AADR v62.0 1240k BIM used to define target positions;
- a PLINK 1240k scaffold with chromosome-position variant identifiers and A/C/G/T
  alleles;
- a two-column chromosome-position-to-rsID map;
- the hs37d5 reference FASTA, 1240k BED and AADR 1240k SNP file;
- the 985-person HGDP/SGDP reference PLINK dataset and its individual-to-population
  mapping, obtained from the original provider.

These controlled or licensed inputs are not redistributed. The public repository
contains only the 48 query coordinates, population-level reference aggregates,
eigenvalues, query manifest and derived metrics needed to verify the reported
calculations and figures. It does not contain raw BAMs, full VCFs, the AADR
individual-level genotype matrix or individual-level reference coordinates.

`pileupCaller --randomHaploid` performs stochastic allele sampling. The recovered
command is exact, but its historical random state was not recorded. The released
processed coordinates and validation files preserve the analysis actually used.

## Query preparation

Use a work directory outside the Git repository because the intermediate PLINK
files are large.

```bash
export PLINK=plink
export BCFTOOLS=bcftools
bash workflows/10_pca/recovery/prepare_imputed_queries.sh \
  query_inputs.tsv AADR_1240K_BIM SCAFFOLD_PREFIX POSITION_TO_RSID \
  OUTPUT_PREFIX WORK_DIR
```

The script preserves the recovered two-stage merge: 16 main treatment datasets
are harmonised first, the eight additional non-UDG datasets are harmonised
separately, and the two groups are then merged to make the 24-query imputed PLINK
dataset.

```bash
export SAMTOOLS=samtools
export PILEUPCALLER=pileupCaller
bash workflows/10_pca/recovery/prepare_pseudohaploid_queries.sh \
  query_inputs.tsv HS37D5_FASTA AADR_1240K_BED AADR_1240K_SNP \
  OUTPUT_PREFIX WORK_DIR
```

The recovered pseudo-haploid command uses `samtools mpileup -B -q30 -Q30 -R`
and `pileupCaller --randomHaploid` across the same 24 treatment datasets.

## Matched-48 projection

After copying the two query PLINK datasets and the reference inputs to the
projection environment, configure the executables in `workflows/config.sh` and
run:

```bash
export ORT_CONFIG="$PWD/workflows/config.sh"
bash workflows/10_pca/recovery/run_recovered_matched48_projection.sh \
  IMPUTED_PREFIX PSEUDOHAPLOID_PREFIX REFERENCE_PREFIX \
  IID_TO_POPULATION workflows/10_pca/recovery/keep_pops_ref.txt OUTPUT_DIR
```

The wrapper sets the genetic-distance column of each query BIM to zero before
PLINK 2 processing, matching the recovered final run and leaving the input PLINK
files unchanged. It then delegates the common harmonisation and projection stages
to the shared workflow and writes the recovered query manifest and build summary.

The final settings were `lsqproject: YES`, `numoutlieriter: 0`,
`numoutevec: 10`, `altnormstyle: NO` and 16 threads. No LD pruning was applied.
The aligned intersection contained 1,134,702 markers; `smartpca` used 1,113,348
after internal filtering. The output contained 985 reference individuals and 48
projected query points, comprising 24 imputed and 24 pseudo-haploid datasets.

## Public validation and figures

Run the canonical checks from the repository root:

```bash
bash workflows/10_pca/recovery/run_validation_and_figures.sh BUILD_DIR
```

The wrapper:

1. recalculates the method-distance table with an absolute tolerance of `1e-12`;
2. validates the reference, query, marker and titration counts;
3. rebuilds Figure 7 and Supplementary Figures S7 and S8 from the released
   query-only coordinates and disclosure-controlled population aggregates;
4. renders each PDF at 150 dpi and, by default, compares its SHA-256 checksum
   with the archived manuscript raster.

The recovery validation passed **20/20 checks**. All three regenerated rasters
were byte-identical to the manuscript rasters. Detailed results are in
`data/processed/pca/validation_results.tsv`, and the publication-safe file map is
`data/processed/pca/input_output_manifest.tsv`.

### Canonical-strict and CI-portable modes

`canonical-strict` is the default and remains unchanged. It requires all 17
numerical and structural checks to pass and all three generated rasters to
match the archived SHA-256 values exactly. The canonical rasters were produced
with Poppler 25.06.0 and TeX Live 2024.

The failing v1.1.0 GitHub Actions runs used Poppler 24.02.0, whose anti-aliasing
produced different PNG bytes even though the underlying numerical outputs and
plotted scientific content were identical. CI therefore selects the portable
mode explicitly:

```bash
ORT_PCA_VALIDATION_MODE=ci-portable \
  bash workflows/10_pca/recovery/run_validation_and_figures.sh BUILD_DIR
```

The portable mode does not weaken numerical validation. It runs the same 17
numerical and structural checks, verifies each bundled canonical raster against
its archived SHA-256 value, requires generated files to be RGB PNGs with the
exact expected dimensions, and compares them with the canonical rasters using
all of these limits:

- normalized mean absolute pixel error no greater than `0.020`;
- fourfold-downsampled grayscale pixel correlation at least `0.995`;
- 32 × 32 difference-hash disagreement no greater than `0.060`;
- bidirectional dark-pixel coverage at least `0.985`, with a 2-pixel tolerance
  for renderer-dependent edge placement.

A missing raster, wrong format, wrong colour mode, wrong dimensions, corrupted
canonical fixture, displaced content, or image outside any threshold fails.
Portable mode tests scientific and visual equivalence without claiming
renderer-dependent byte identity; it does not replace the canonical strict
record.
