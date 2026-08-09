# 08 Kinship estimation

The study uses complementary kinship summaries rather than treating them as fully
independent analyses.

## READv2

The primary reported result came from an archived ORT1--ORT27 1240k callset that
included chromosomes 1--22, X, and Y. `run_readv2_primary_chrxy.sh` reconstructs
the recorded calling logic with `pileupCaller --randomHaploid` and READv2's
within-cohort median baseline. READv2 retained 21 individuals in the archived
output. `run_readv2_cemetery.sh` reconstructs a later ORT6--ORT27 autosome-only
sensitivity run, which returned the same first-degree parent--offspring subtype.

For both runs, each cemetery input combined full-UDG libraries a--c and was
soft-clipped by 6 bp at both read ends before pseudo-haploid calling. The public
wrappers accept these already prepared BAMs through manifests; they do not recreate
the upstream library merge or terminal soft-clipping. The complete cemetery cohorts
are not available through PRJEB112497 alone.

```bash
ORT_CONFIG=/path/to/config.sh bash workflows/08_kinship/run_readv2_primary_chrxy.sh \
  workflows/08_kinship/readv2_primary_bam_manifest.example.tsv output/readv2_primary

ORT_CONFIG=/path/to/config.sh bash workflows/08_kinship/run_readv2_cemetery.sh \
  workflows/08_kinship/readv2_bam_manifest.example.tsv output/readv2_autosomes
```

## KING and IBS0

`run_king_ibs0.sh` merges the two full-UDG GLIMPSE VCFs and restricts sites to
reference-panel allele frequencies between 0.1 and 0.9.
`summarize_king_ibs0.py` calculates KING-robust kinship, opposing homozygotes,
the unrelated expectation `mean(2 p^2 q^2)`, and the full-sibling expectation
(one quarter of the unrelated value). The same calculation can be restricted to
high-posterior-probability genotypes.

The analysis requires individual-level BAM/VCF inputs and the relevant reference
frequency panel; those files are not distributed.
