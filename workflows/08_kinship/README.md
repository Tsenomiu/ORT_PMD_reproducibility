# 08 Kinship estimation

The study uses complementary kinship summaries rather than treating them as fully
independent analyses.

## READv2

`run_readv2_cemetery.sh` creates pseudo-haploid calls at 1240k sites with
`pileupCaller --randomHaploid`, restricts to autosomes, and runs READv2 using the
within-cohort median baseline. The complete 22-person normalisation cohort is not
available through PRJEB112497 alone.

## KING and IBS0

`run_king_ibs0.sh` merges the two full-UDG GLIMPSE VCFs and restricts sites to
reference-panel allele frequencies between 0.1 and 0.9.
`summarize_king_ibs0.py` calculates KING-robust kinship, opposing homozygotes,
the unrelated expectation `mean(2 p^2 q^2)`, and the full-sibling expectation
(one quarter of the unrelated value). The same calculation can be restricted to
high-posterior-probability genotypes.

The analysis requires individual-level BAM/VCF inputs and the relevant reference
frequency panel; those files are not distributed.
