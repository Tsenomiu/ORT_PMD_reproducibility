# Data access

## Sequencing reads

Study reads are deposited in the European Nucleotide Archive under accession
**PRJEB112497**. ENA is the authoritative source for filenames, checksums, metadata,
and access status. The study remains private during preparation of this local
release candidate; the accession should be described as publicly accessible only
after its ENA release has been verified.

| Library | ENA sample | BioSample | Library type |
|---|---|---|---|
| ORT15_a | ERS30022076 | SAMEA122427376 | full-UDG |
| ORT15_b | ERS30022077 | SAMEA122427377 | full-UDG |
| ORT15_c | ERS30022078 | SAMEA122427378 | full-UDG |
| ORT15_nu | ERS30022079 | SAMEA122427379 | non-UDG |
| ORT16_a | ERS30022080 | SAMEA122427380 | full-UDG |
| ORT16_b | ERS30022081 | SAMEA122427381 | full-UDG |
| ORT16_c | ERS30022082 | SAMEA122427382 | full-UDG |
| ORT16_nu | ERS30022083 | SAMEA122427383 | non-UDG |

This repository does not mirror FASTQ files.

## Included data

The repository contains only compact inputs needed to rebuild figures or verify
reported calculations, including:

- per-bin imputation summaries and chromosome-jackknife counts;
- ORT15 and ORT16 PCA coordinates, eigenvalues, and population-level reference summaries;
- filtered ORT15–ORT16 ancIBD plotting tables;
- mapDamage text summaries used for PMD plots; and
- the complete 16-state TKGWV2 result table, the complete 210-pair primary READv2
  output, anonymized READv2 normalization vectors, ORT15–ORT16 sensitivity rows,
  and compact mitochondrial summaries;
- publication-safe PCA recovery manifests and Table 1 generation/provenance
  records; and
- small tabular inputs used by the kinship and imputation figures.

These files are derived results, not raw sequence data.

## Data not included

- BAM/CRAM, VCF/BCF, and HDF5 files;
- PLINK or EIGENSTRAT genotype matrices;
- hs37d5, 1000 Genomes, AADR, HapMapChrX, or schmutzi reference resources;
- exact archaeological coordinates or the original field-report site plan;
- the underlying cemetery-wide READv2 BAM cohort and the full autosome-only cohort
  output; and
- raw micro-CT volumes.

## Analysis-specific access boundaries

| Analysis | Data route | Additional requirement |
|---|---|---|
| Read processing, PMD correction, imputation | ENA PRJEB112497 | Reference genome, panels, and server-scale compute |
| Concordance and jackknife | Rebuild VCFs from ENA; compact counts included | 1000 Genomes frequencies and GLIMPSE resources |
| ancIBD | Compact ORT15–ORT16 plotting inputs included | AADR v62 and 1000 Genomes resources for a full rerun |
| READv2 | Complete primary 210-pair output, anonymized normalization vectors, and focal sensitivity rows included | Controlled BAMs and the unrecovered primary pseudo-haploid random state prevent byte-identical raw-input reconstruction |
| PCA | ORT15 and ORT16 coordinates and aggregate backgrounds included | AADR genotype data for exact reference reconstruction |
| Site map | Figure available with the article | Controlled archaeological spatial information |
| Micro-CT | Export code and acquisition metadata included | Original VOX volumes and selected slice coordinates |

Requests for controlled archaeological or imaging data should be directed to the
corresponding author and the relevant curating institution. Access may require
institutional approval.
