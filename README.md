# ORT post-mortem damage correction workflows

Analysis workflows, figure generators, and compact derived data for:

> **Trade-offs in post-mortem damage correction for low-coverage ancient genomes:
> site retention, imputation and kinship inference**

The study compares terminal trimming, base-quality rescaling, and SNP masking in
matched full-UDG and non-UDG libraries from two ancient individuals, ORT15 and
ORT16.

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
| `figures/` | Generators for main Figures 1–7 and Supplementary Figures S1–S8 |
| `data/summary/` | Compact derived tables used by figures and validation checks |
| `docs/` | Detailed methods, data-access boundaries, external resources, and limitations |
| `software_versions.tsv` | Recorded software versions and their roles |

`reproduce_all.sh` checks the included summaries, reconstructs the paired
chromosome jackknife and PCA metrics, validates the TKGWV2 and READv2 tables,
and rebuilds the supported figures.

| Component | What runs from included files |
|---|---|
| Numerical checks | Concordance, jackknife, TKGWV2, READv2, mitochondrial, and PCA summaries |
| Figures | Main Figures 1–7 and Supplementary Figures S1 and S4–S8 |
| External-data workflows | Portable commands are supplied; raw or licensed inputs are required |
| Supplementary Figures S2–S3 | Code and metadata are supplied; controlled source data are required |
| TKGWV2 | The complete 16-state result table is validated; no unverified runner is supplied |

## What requires external data

The complete workflows are provided, but an end-to-end rerun from raw reads
requires separately obtained sequencing data, reference resources,
configuration, and suitable compute. The repository does not contain FASTQ,
BAM, VCF, HDF5, third-party reference panels, exact archaeological coordinates,
or raw micro-CT volumes.

Raw sequencing reads are deposited in the European Nucleotide Archive under
study accession **PRJEB112497**. Third-party resources such as hs37d5, 1000
Genomes Phase 3, AADR v62, HapMapChrX, and the schmutzi mitochondrial panel must
be obtained from their original providers.

Supplementary Figure S2 requires controlled archaeological spatial inputs.
Supplementary Figure S3 requires the original CT volumes and recorded slice
coordinates. The PCA figures included here use population-level reference
aggregates; reconstructing the exact individual-level AADR background requires
separately licensed AADR data.

See [Data access](docs/DATA_ACCESS.md),
[Third-party resources](docs/THIRD_PARTY_RESOURCES.md), and
[Known limitations](docs/KNOWN_LIMITATIONS.md) before attempting a full rerun.

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

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Please also
cite the associated article and the primary software and resource papers listed
in [`docs/THIRD_PARTY_RESOURCES.md`](docs/THIRD_PARTY_RESOURCES.md).

## License

This repository is released under the [MIT License](LICENSE).
