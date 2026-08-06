# ORT post-mortem damage correction workflows

Analysis and figure code for the study:

> **Post-mortem damage correction in low-coverage ancient genomes: effects on
> site retention, imputation and kinship inference**

The study compares terminal trimming, base-quality rescaling, and SNP masking in
matched full-UDG and non-UDG libraries from two ancient individuals (ORT15 and
ORT16). This repository contains portable workflows, deterministic plotting code,
and the compact inputs required for the included numerical checks.

## Repository contents

| Path | Description |
|---|---|
| [`workflows/`](workflows/) | Numbered analysis stages from read processing through PCA |
| [`figures/`](figures/) | Plotting scripts and minimal figure inputs |
| [`data/summary/`](data/summary/) | Compact tables used by figures and numerical tests |
| [`docs/METHODS_DETAILS.md`](docs/METHODS_DETAILS.md) | Additional computational-method details |
| [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md) | Sequencing-data access and exclusions |
| [`docs/THIRD_PARTY_RESOURCES.md`](docs/THIRD_PARTY_RESOURCES.md) | External reference resources required by the workflows |
| [`docs/KNOWN_LIMITATIONS.md`](docs/KNOWN_LIMITATIONS.md) | Analyses that require controlled or third-party inputs |
| [`software_versions.tsv`](software_versions.tsv) | Software versions used in the study |

The repository does not contain raw sequencing data, BAM/VCF intermediates,
reference panels, exact archaeological coordinates, or raw CT volumes.

## Reproduce the included results

Python 3.11 or later is recommended.

```bash
git clone https://github.com/Tsenomiu/ORT_PMD_reproducibility.git
cd ORT_PMD_reproducibility
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
bash check_repository.sh
bash reproduce_all.sh
```

Generated files are written under `_build/`, which is ignored by Git. The local
reproduction command:

- checks shell and Python syntax;
- reconstructs the paired chromosome-jackknife results;
- verifies the included PCA summary metrics; and
- rebuilds figures whose compact inputs can be distributed.

The site map requires controlled archaeological spatial data. The micro-CT exporter
requires the original VOX volumes and selected slice coordinates. Exact PCA reference
backgrounds require separately licensed AADR inputs. These external-input boundaries
are described in [`docs/KNOWN_LIMITATIONS.md`](docs/KNOWN_LIMITATIONS.md).

## Re-run the analysis workflows

The numbered workflows require substantial compute and data obtained separately.
Create an untracked local configuration:

```bash
cp workflows/config.example.sh workflows/config.sh
# Edit paths to local software, reads, and reference resources.
export ORT_CONFIG="$PWD/workflows/config.sh"
```

Then follow the README in each numbered directory. Example manifests describe the
required columns without embedding local server paths. Scripts that need BAM, VCF,
reference-genome, 1000 Genomes, AADR, or schmutzi resources do not download or
redistribute those files.

## Data access

Sequencing reads are deposited in the European Nucleotide Archive under study
accession **PRJEB112497**. Obtain authoritative filenames and checksums directly from
ENA. Reference resources such as hs37d5, 1000 Genomes Phase 3, AADR v62, HapMapChrX,
and the schmutzi mitochondrial panel must be obtained from their original providers.

See [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md) and
[`docs/THIRD_PARTY_RESOURCES.md`](docs/THIRD_PARTY_RESOURCES.md).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Cite the associated
article and the versioned Zenodo archive when they become available. GitHub contains
the maintained code; Zenodo will provide immutable snapshots of reviewed releases.

## License

No public reuse license has yet been selected. The repository will remain private
until the authors add explicit licenses for code and distributable derived data.
