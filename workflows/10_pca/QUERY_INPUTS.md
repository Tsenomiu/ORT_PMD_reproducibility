# PCA query inputs

The projection workflow expects three PLINK datasets:

- 24 GLIMPSE-imputed ORT treatment datasets;
- the same 24 datasets called pseudo-haploid at AADR 1240k sites with
  `pileupCaller --randomHaploid`; and
- 985 present-day HGDP/SGDP AADR-v62 reference individuals from 154 populations.

The 24 treatments comprise four full-UDG and eight non-UDG treatments per
individual. Build one PLINK trio per query representation, confirm sample order
against `data/summary/pca/query_manifest_matched48.tsv`, and supply both trios to
`run_matched48_projection.sh` with the reference trio.

This repository begins the projection workflow from those three prebuilt PLINK
datasets. It does not include a verified FASTQ/BAM-to-query-trio driver for all 24
imputed and 24 pseudo-haploid datasets. The recorded pseudo-haploid calls used
pileupCaller (sequenceTools) 1.5.4.0; the missing component is the complete upstream
driver, not the software version. This boundary is also listed in
`docs/KNOWN_LIMITATIONS.md`.

AADR genotypes and the individual-level reference coordinates are obtained from the
provider and are not stored in this repository.
