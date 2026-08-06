# Known reproducibility limitations

This repository supports deterministic regeneration of included summary statistics
and figures from compact inputs. It is not a self-contained copy of the full analysis
environment.

## External data and compute

- Read processing, PMD correction, imputation, and most kinship workflows require
  study reads from ENA accession PRJEB112497 and separately obtained reference data.
- BAM, VCF/BCF, HDF5, PLINK, and EIGENSTRAT intermediates are not distributed.
- The complete 22-person READv2 normalisation cohort is not contained in
  PRJEB112497, so that cohort-level normalisation cannot be rerun from this repository
  alone.
- AADR and 1000 Genomes resources retain their provider terms and are not bundled.

## Figures requiring non-distributed inputs

- The archaeological site map requires controlled spatial data and the authorised
  site-plan source.
- The micro-CT exporter requires the original VOX volumes and selected slice
  coordinates.
- PCA plotting code includes aggregate reference backgrounds suitable for checking
  the focal ORT coordinates. Reproducing the exact individual-level reference
  background requires separately licensed AADR data.

## Workflow coverage

The numbered workflow directories provide final scripts or portable command wrappers
for the main analysis stages. Each README lists the required inputs and expected
outputs. A TKGWV2 runner is not included because the exact locally patched conversion
code used in the study is not available as a verified, portable script; the repository
does not substitute an untested reconstruction.

Local automated checks therefore test syntax, compact numerical reconstructions, and
figure generation. They do not execute server-scale processing from FASTQ files.
