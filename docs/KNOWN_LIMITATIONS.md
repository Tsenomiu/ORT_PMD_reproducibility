# Known reproducibility limitations

This repository supports deterministic regeneration of included summary statistics
and figures from compact inputs. It is not a self-contained copy of the full analysis
environment.

## External data and compute

- Read processing, PMD correction, imputation, and most kinship workflows require
  study reads from ENA accession PRJEB112497 and separately obtained reference data.
- BAM, VCF/BCF, HDF5, PLINK, and EIGENSTRAT intermediates are not distributed.
- The complete 210-pair output from the primary READv2 cohort is included as a
  derived table. The input BAM cohort and the later autosome-only cohort inputs are
  controlled and are not contained in PRJEB112497, so cohort-level normalisation
  cannot be rerun from this repository alone.
- AADR and 1000 Genomes resources retain their provider terms and are not bundled.

## Figures requiring non-distributed inputs

- The archaeological site map requires controlled spatial data and the authorised
  site-plan source.
- The micro-CT exporter requires the original VOX volumes and selected slice
  coordinates.
- PCA plotting code includes aggregate reference backgrounds suitable for checking
  the ORT15 and ORT16 coordinates. Reproducing the exact individual-level reference
  background requires separately licensed AADR data.
- The matched PCA projection wrapper starts from prebuilt PLINK query datasets. The
  repository does not reconstruct all 24 imputed and 24 pseudo-haploid query datasets
  from FASTQ/BAM inputs, and the exact pileupCaller build used for pseudo-haploid
  calling was not recorded.

## Workflow coverage

The numbered workflow directories provide final scripts or portable command wrappers
for the main analysis stages. Each README lists the required inputs and expected
outputs. The complete 16-state TKGWV2 result table is included and validated, but a
TKGWV2 runner is not supplied because the exact locally patched conversion code used
in the study is not available as a verified, portable script; the repository does not
substitute an untested reconstruction.

Local automated checks therefore test syntax, compact numerical reconstructions, and
figure generation. They do not execute server-scale processing from FASTQ files.
