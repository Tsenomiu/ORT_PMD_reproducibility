# ORT PMD reproducibility workflow

Analysis workflows, data-figure generators, and compact derived data for:

> **Matched full-UDG and non-UDG ancient DNA libraries reveal trade-offs in
> post-mortem damage correction for imputation and kinship inference**

The study compares terminal trimming, base-quality rescaling, and SNP masking in
matched full-UDG and non-UDG libraries from two ancient individuals, ORT15 and
ORT16. Release **v1.2.0** corrects the mitochondrial input workflow by using one
collapsed-only full-UDG read-set BAM per individual and prohibiting overlapping
nested-state merges. D4o1 assignments are unchanged; corrected caller-emitted SNP
counts are 36 for ORT15 and 37 for ORT16, with 36 exact shared calls. The release
also adds the Reduce18F/Round45 DeltaALT genomic-block uncertainty summaries and
confidence intervals used by Figure 2.

## Quick start

The commands below verify the included numerical summaries and rebuild every
figure supported by the compact data in this repository.

Run them from the repository root.

```bash
git clone https://github.com/Tsenomiu/ORT_PMD_reproducibility.git
cd ORT_PMD_reproducibility
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
bash reproduce_all.sh
```

Results are written to `_build/`. A successful run ends with:

```text
Reproduction complete: all locally reproducible checks passed.
```

The default PCA check is `canonical-strict`: it requires the three regenerated
150-dpi rasters to match the archived SHA-256 values exactly. Those canonical
rasters were produced with Poppler 25.06.0 and TeX Live 2024. GitHub Actions
uses the explicitly selected `ci-portable` mode, which keeps every numerical,
coordinate, source-data, and plotting-input check strict while checking raster
format, dimensions, and documented visual-equivalence thresholds. See the
[PCA recovery workflow](workflows/10_pca/recovery/README.md) for details.

`pdfTeX` (from TeX Live or MacTeX) is needed to create the final PDF wrappers for
Supplementary Figures S7 and S8. If it is unavailable, the script completes the
remaining checks and reports those two figures as skipped.

Run only the repository safety and syntax checks with:

```bash
bash check_repository.sh
```

## What is included

| Path | Contents |
|---|---|
| `workflows/` | Ten analysis stages, from sequence processing through PCA |
| `figures/` | Generators for main Figures 2–7 and Supplementary Figures S1–S8 |
| `data/summary/` | Compact derived tables used by figures and validation checks |
| `data/processed/` | Curated processed verification data for recovered workflows |
| `docs/` | Detailed methods, data-access boundaries, external resources, and limitations |
| `software_versions.tsv` | Recorded software versions and their roles |

`reproduce_all.sh` checks the included summaries, reconstructs the paired
chromosome jackknife and PCA metrics, validates the TKGWV2 and READv2 tables,
and rebuilds the supported figures.

| Component | What runs from included files |
|---|---|
| Numerical checks | DeltaALT uncertainty, concordance, jackknife, TKGWV2, READv2, corrected mitochondrial, and PCA summaries |
| Figures | Main Figures 2–7 and Supplementary Figures S1 and S4–S8 |
| External-data workflows | Portable commands are supplied; raw or licensed inputs are required |
| Supplementary Figures S2–S3 | Code and metadata are supplied; controlled source data are required |
| TKGWV2 | Recovered patch, processed-input runner, and all 16 exact validation rows |
| READv2 | Derived normalization inputs and validation of the reported primary and sensitivity results |
| PCA | Path-neutral recovery workflow, query-only coordinates, aggregate references, canonical exact raster checks, and CI-portable visual checks |
| Table 1 | Original generator, output, and 80-cell provenance validation |

Figure 1 is an author-prepared workflow schematic created in PowerPoint, not a
data-derived figure. Its source artwork is therefore not included in this
repository and `reproduce_all.sh` does not regenerate it.

## What requires external data

The complete workflows are provided, but an end-to-end rerun from raw reads
requires separately obtained sequencing data, reference resources,
configuration, and suitable compute. The repository does not contain FASTQ,
BAM, VCF, HDF5, third-party reference panels, exact archaeological coordinates,
or raw micro-CT volumes.

Raw sequencing reads are publicly available in the European Nucleotide Archive
under study accession **PRJEB112497**. ENA is the authoritative source for
filenames, checksums, metadata, and current access status. The
files in `data/summary/` and `data/processed/` are compact processed or derived
verification data, not raw sequence data.

Third-party resources such as hs37d5, 1000
Genomes Phase 3, AADR v62, HapMapChrX, and the schmutzi mitochondrial panel must
be obtained from their original providers.

Supplementary Figure S2 requires controlled archaeological spatial inputs.
Supplementary Figure S3 requires the original CT volumes and recorded slice
coordinates. The PCA figures included here use population-level reference
aggregates; reconstructing the exact individual-level AADR background requires
separately licensed AADR data.

See [Data access](docs/DATA_ACCESS.md),
[Third-party resources](docs/THIRD_PARTY_RESOURCES.md), and
[Reproducibility limitations](docs/reproducibility_limitations.md) before
attempting a full rerun.

## Recovered workflow entry points

- TKGWV2: `workflows/08_kinship/tkgwv2/`
- READv2: `workflows/08_kinship/readv2/`
- Corrected mitochondrial calling: `workflows/09_mtdna/`
- DeltaALT uncertainty: `workflows/04_alt_fraction/uncertainty/`
- PCA recovery: `workflows/10_pca/recovery/`
- Table 1: `workflows/06_concordance/table1/`
- Main-output source map: `docs/main_output_source_manifest.tsv`
- v1.1.0 scientific validation record: `v1.1.0_VALIDATION_REPORT.md`
- Release exclusions: `v1.1.0_EXCLUSIONS.md`
- v1.1.1 maintenance notes: `v1.1.1_RELEASE_NOTES.md`
- v1.2.0 correction and uncertainty notes: `v1.2.0_RELEASE_NOTES.md`
- v1.2.0 correction validation: `v1.2.0_VALIDATION_REPORT.md`
- v1.2.0 changed-file inventory: `v1.2.0_FILE_MANIFEST.tsv`

The exact original primary BAM-to-genotype command and random state used for
pseudo-haploid sampling in READv2 were not preserved. The derived genotype
inputs, downstream READv2 commands, normalization procedure, and reported
outputs are provided to permit verification of the published calculations.

## Workflow configuration

Copy the example configuration and replace every placeholder with a local path:

```bash
cp workflows/config.example.sh workflows/config.sh
# Edit workflows/config.sh, then load its variables:
export ORT_CONFIG="$PWD/workflows/config.sh"
source "$ORT_CONFIG"
```

Then follow the input and output instructions in the relevant workflow README.
Executable examples are supplied where the complete command was reconstructed.
Large genomic inputs and outputs are ignored by `.gitignore` and should remain
outside the repository.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). When citing
the associated manuscript, use the title **Matched full-UDG and non-UDG ancient
DNA libraries reveal trade-offs in post-mortem damage correction for imputation
and kinship inference**. Please also cite the primary software and resource
papers listed in
[`docs/THIRD_PARTY_RESOURCES.md`](docs/THIRD_PARTY_RESOURCES.md).

## License

Author-written code and derived tables are released under the
[MIT License](LICENSE). The TKGWV2 compatibility patch is a modification of
GPL-2.0 software and is distributed under GPL-2.0-only; see
[`workflows/08_kinship/tkgwv2/THIRD_PARTY_NOTICE.md`](workflows/08_kinship/tkgwv2/THIRD_PARTY_NOTICE.md).
