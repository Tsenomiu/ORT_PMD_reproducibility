# ORT PMD reproducibility workflow

Analysis scripts, figure generators and compact derived data for:

> **Matched full-UDG and non-UDG ancient DNA libraries reveal trade-offs in
> post-mortem damage correction for imputation and kinship inference**

The study compares trimming, base-quality rescaling and SNP masking in matched
libraries from ORT15 and ORT16.

## Quick start

```bash
git clone https://github.com/Tsenomiu/ORT_PMD_reproducibility.git
cd ORT_PMD_reproducibility
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
ORT_PCA_VALIDATION_MODE=ci-portable bash reproduce_all.sh
```

This verifies the included numerical summaries and rebuilds supported figures
in `_build/`. PCA raster validation and the titration PDFs also require
`pdftex` and `pdftoppm`; unavailable renderers produce explicit skips.

`ci-portable` uses the same numerical checks as the default `canonical-strict`
mode but allows documented renderer-dependent raster differences. For exact
canonical raster checks, run `bash reproduce_all.sh` with the recorded rendering
environment. See [validation](docs/VALIDATION.md) for details.

## Contents

| Directory | Purpose |
|---|---|
| [workflows](workflows/README.md) | Ten analysis stages, configuration and input instructions |
| [figures](figures/README.md) | Data-figure generators and input requirements |
| [data/summary](data/summary/README.md) | Derived tables used by figures and numerical checks |
| [data/processed](data/processed/README.md) | Processed verification inputs and records |
| [docs](docs/METHODS_DETAILS.md) | Methods, validation and reproducibility boundaries |

Figure 1 was created in PowerPoint and is not regenerated here. The compact
package rebuilds main Figures 2–7 and Supplementary Figures S1–S4 and S7–S8;
S5–S6 require external spatial or CT inputs. See the [figure map](figures/README.md).

## Data and reproducibility

Raw reads are public in ENA under
[PRJEB112497](https://www.ebi.ac.uk/ena/browser/view/PRJEB112497).
BAMs, VCFs, reference panels, exact archaeological coordinates and raw CT volumes
are not bundled. A raw-read rerun requires external inputs and suitable compute.
The compact checks do not claim byte-identical regeneration of the historical
pseudo-haploid calls.

- [Data access](docs/DATA_ACCESS.md) and [third-party resources](docs/THIRD_PARTY_RESOURCES.md)
- [Reproducibility limitations](docs/reproducibility_limitations.md)
- [Validation and workflow provenance](docs/VALIDATION.md)
- [Software versions](software_versions.tsv) and [changelog](CHANGELOG.md)

## Citation and licence

Use [CITATION.cff](CITATION.cff) and cite the relevant software/resource papers.
Author-written code and derived tables use the [MIT licence](LICENSE).
The TKGWV2 compatibility patch is GPL-2.0-only; see its
[third-party notice](workflows/08_kinship/tkgwv2/THIRD_PARTY_NOTICE.md).
