# READv2 verification and rerun package

This directory documents the three accepted READv2 analyses used to assess the
ORT15–ORT16 relationship. It separates exact recovery from reconstruction and
does not claim complete primary raw-BAM-to-result reproducibility.

> **Important limitation:** The exact original primary BAM-to-genotype command and random state used for pseudo-haploid sampling were not preserved. The derived genotype inputs, downstream READv2 commands, normalization procedure and reported outputs are provided to permit verification of the published calculations.

For this public package, the safe verification inputs are derived, non-genomic
P0 tables; full PLINK genotype matrices are not distributed. The primary
archived genotype input is retained in the author archive. The path-neutral
primary wrapper starts from its preserved 27-person PLINK BED/BIM/FAM boundary.
The later autosome-only and transversion-only wrapper records the available
genotype-generation commands, but it requires controlled cemetery BAMs and
third-party reference files. Its recovered command used `--randomHaploid`
without an explicit seed, so a new BAM-level run need not be byte-identical to
the archived one.

## Software identity

READv2 is from [GuntherLab/READv2](https://github.com/GuntherLab/READv2), tag
`2.1.0`, commit `73739bef65b3725315417e70ad8edf6aca21f7b3`. The executed
`READ2.py` prints version `v2.01` and has SHA-256
`63e8915a9887b756173547d0655611a796450cc5cef28c17f79daff4ad1e3134`.
The upstream program is GPL-3.0 licensed and is not vendored in this MIT-licensed
repository. Obtain that exact revision from the upstream project.

Recorded software environments are listed in `environment.tsv`.

## Analysis definitions and exact commands

### Primary chromosomes 1–22, X and Y

The archived primary PLINK input contained ORT1–ORT27. These are the exact
recovered downstream commands, with paths made configurable by
`run_primary_from_processed_plink.sh`:

```bash
plink --bfile ORT_fix --mind 0.8 --make-bed --out filtered
python READ2.py -i filtered --window_size 500000
```

The `--mind 0.8` filter removed ORT1, ORT17, ORT19, ORT20, ORT22 and ORT25.
Twenty-one people remained, giving 210 pairwise comparisons. The explicit
500-kb window option is part of the original command.

To rerun from the preserved processed input boundary:

```bash
READ2_PY=/path/to/READv2/READ2.py \
  bash workflows/08_kinship/readv2/run_primary_from_processed_plink.sh \
  /path/to/ORT_fix output/readv2_primary
```

### Autosome-only 1240k and transversion-only sensitivities

The recovered later workflow jointly processed ORT6–ORT27. It used one pileup,
created an all-1240k pseudo-haploid callset and a matching transversion-only
callset, restricted both to chromosomes 1–22, and ran READv2 with its default
5-Mb windows. Here “all-1240k” means all available targets on the
AADR v62.0 1240k panel, including transitions; it does not mean every site in a
whole-genome reference panel.

The scientifically relevant options are preserved verbatim in
`run_autosome_and_transversion_sensitivities.sh`:

```text
samtools mpileup -B -q30 -Q30 -R -l PANEL.bed -f hg19.fa <22 BAMs>
pileupCaller --randomHaploid --sampleNames ORT6,...,ORT27 --samplePopName ORT -f PANEL.snp -p cemetery_all
pileupCaller --randomHaploid --skipTransitions --sampleNames ORT6,...,ORT27 --samplePopName ORT -f PANEL.snp -p cemetery_tv
plink --bfile cemetery_all --chr 1-22 --make-bed --allow-no-sex --threads 8 --memory 20000 --out read2_allsites/cemetery
python3 READ2.py -i read2_allsites/cemetery
plink --bfile cemetery_tv --chr 1-22 --make-bed --allow-no-sex --threads 8 --memory 20000 --out read2_transversions/cemetery
python3 READ2.py -i read2_transversions/cemetery
```

Prepare a two-column manifest from `bam_manifest.example.tsv`, then run:

```bash
HG19_FASTA=/path/to/hg19.fa \
AADR_1240K_BED=/path/to/1240k.all.chr.bed \
AADR_V62_SNP=/path/to/v62.0_1240k_public.snp \
READ2_PY=/path/to/READv2/READ2.py \
  bash workflows/08_kinship/readv2/run_autosome_and_transversion_sensitivities.sh \
  /path/to/bam_manifest.tsv output/readv2_sensitivities
```

The BAMs used in the original later run were merged full-UDG libraries after
6-bp terminal soft clipping. Those controlled inputs and third-party reference
files are not distributed in this repository.

## Normalization and reported results

READv2 used its default within-cohort median normalization:

```text
baseline = median(raw P0 across every cohort pair)
P0_norm  = focal raw P0 / baseline
KC       = 1 - P0_norm
```

The exact archived focal values and the manuscript-level rounded values are:

| Analysis | P0_norm | Rounded P0_norm | KC | Rounded KC | Overlap SNPs | Degree | Subtype |
|---|---:|---:|---:|---:|---:|---|---|
| Primary chr1–22/X/Y | 0.7690333381022743 | 0.769 | 0.23096666189772574 | 0.231 | 169,172 | First degree | Parent–offspring |
| Autosome-only 1240k | 0.7684078605188551 | **0.768** | 0.23159213948114488 | **0.232** | 157,451 | First degree | Parent–offspring |
| Transversion-only | 0.7699815124598774 | 0.770 | 0.23001848754012255 | **0.230** | 31,100 | First degree | N/A |

The transversion-only result supports first-degree relatedness but does not meet
READv2's within-first-degree parent–offspring rule. It must not be described as
an additional parent–offspring subtype call.

## Included publication files

| File | Classification | Purpose |
|---|---|---|
| `run_primary_from_processed_plink.sh` | publication wrapper around original recovered commands | Exact primary downstream filtering and READv2 options from the preserved PLINK boundary |
| `run_autosome_and_transversion_sensitivities.sh` | path-neutral publication copy of original recovered code | Available later genotype generation and exact READv2 commands |
| `bam_manifest.example.tsv` | publication template | Required ordering and columns for controlled later-run BAM inputs |
| `environment.tsv` | recovered environment metadata | Exact recorded versions and roles |
| `validate_reported_results.py` | generated validation code | Recomputes normalization, KC, calls, sample counts and required rounded values |
| `data/processed/readv2/focal_results.tsv` | processed source data | Exact focal output rows for all three accepted analyses |
| `data/processed/readv2/*_p0_anonymized.tsv` | processed source data | Pairwise P0 values with unrelated pair identities replaced by numerical ranks |
| `data/processed/readv2/sample_manifest.tsv` | processed source data | Complete analysis membership and the six primary `--mind` exclusions |
| `data/processed/readv2/analysis_summary.tsv` | generated summary | Exact and three-decimal focal results |
| `data/processed/readv2/validation_results.tsv` | generated validation output | Machine-readable PASS/FAIL record |

No raw BAM, FASTQ, VCF or full PLINK genotype matrix is included. The
anonymized P0 tables retain every value needed to recompute the three cohort
medians while avoiding redistribution of unrelated cemetery pair identities.

## Validate locally

The validator uses only the Python standard library:

```bash
python3 workflows/08_kinship/readv2/validate_reported_results.py
```

It rewrites `data/processed/readv2/validation_results.tsv` deterministically and
exits non-zero on any mismatch.
