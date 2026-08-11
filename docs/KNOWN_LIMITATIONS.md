# Known reproducibility limitations

The release-candidate limitation statement is maintained in
[`reproducibility_limitations.md`](reproducibility_limitations.md). The summary
below is retained for compatibility with links from v1.0.1.

This repository supports deterministic regeneration of included summary statistics
and figures from compact inputs. It is not a self-contained copy of the full analysis
environment.

## External data and compute

- Read processing, PMD correction, imputation, and most kinship workflows require
  study reads from ENA accession PRJEB112497 and separately obtained reference data.
- BAM, VCF/BCF, HDF5, PLINK, and EIGENSTRAT intermediates are not distributed.
- Derived READv2 normalization inputs and focal results support recomputation of
  the reported values. The exact primary BAM-to-genotype command and pseudo-haploid
  random state were not preserved, so raw-BAM-to-result byte identity is not
  claimed.
- AADR and 1000 Genomes resources retain their provider terms and are not bundled.

## Figures requiring non-distributed inputs

- The archaeological site map requires controlled spatial data and the authorised
  site-plan source.
- The micro-CT exporter requires the original VOX volumes and selected slice
  coordinates.
- PCA plotting code includes aggregate reference backgrounds suitable for checking
  the ORT15 and ORT16 coordinates. Reproducing the exact individual-level reference
  background requires separately licensed AADR data.
- The recovered PCA package documents TEST query preparation and the authoritative
  RC projection. Licensed individual-level reference inputs and controlled ORT
  genomic inputs remain external. The recorded pseudo-haploid calls used
  pileupCaller (sequenceTools) 1.5.4.0.

## Workflow coverage

The numbered workflow directories provide final scripts or portable command wrappers
for the main analysis stages. Each README lists the required inputs and expected
outputs. The recovered TKGWV2 package supplies the exact four-line compatibility
patch, a path-neutral runner reconstructed from successful logs, and exact validation
of all 16 reported rows from the retained processed-input boundary.

Local automated checks therefore test syntax, compact numerical reconstructions, and
figure generation. They do not execute server-scale processing from FASTQ files.
