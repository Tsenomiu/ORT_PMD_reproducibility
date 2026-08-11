# 08 Kinship estimation

The study uses complementary kinship summaries rather than treating them as fully
independent analyses.

The proposed v1.1.0 recovery additions are organized under:

- [`tkgwv2/`](tkgwv2/) for the recovered TKGWV2 patch, processed-input runner,
  manifests, and exact 16-row validation;
- [`readv2/`](readv2/) for the recovered downstream commands, path-neutral
  verification workflow, normalization data, and limitation statement.

## READv2

The complete curated documentation and verification entry point is
[`readv2/`](readv2/). It distinguishes the primary chromosomes 1--22/X/Y
analysis from the later autosome-only and transversion-only sensitivities and
provides the cohort-median normalization values needed to recompute the reported
results.

```bash
python3 workflows/08_kinship/readv2/validate_reported_results.py
```

The exact original primary BAM-to-genotype command and random state used for
pseudo-haploid sampling were not preserved. The package therefore verifies the
reported values from recovered derived/downstream records but does not claim
byte-identical reconstruction from the primary BAMs.

## KING and IBS0

`run_king_ibs0.sh` merges the two full-UDG GLIMPSE VCFs and restricts sites to
reference-panel allele frequencies between 0.1 and 0.9.
`summarize_king_ibs0.py` calculates KING-robust kinship, opposing homozygotes,
the unrelated expectation `mean(2 p^2 q^2)`, and the full-sibling expectation
(one quarter of the unrelated value). The wrapper generates both the genome-wide
summary and a chromosome-1 sensitivity summary requiring maximum genotype
probability of at least 0.99 in both individuals.

```bash
ORT_CONFIG=/path/to/config.sh bash workflows/08_kinship/run_king_ibs0.sh \
  /path/to/ORT15.full_udg.imputed.vcf.gz \
  /path/to/ORT16.full_udg.imputed.vcf.gz \
  output/king_ibs0
```

The named outputs are `king_ibs0_genomewide.tsv` and
`king_ibs0_chr1_gp99.tsv`.

The analysis requires individual-level BAM/VCF inputs and the relevant reference
frequency panel; those files are not distributed.
