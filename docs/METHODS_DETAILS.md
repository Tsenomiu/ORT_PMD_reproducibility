# Additional methods

This document provides computational details that complement the associated
article. Executable commands and input schemas are organised by analysis stage in
[`../workflows/`](../workflows/).

## Library preparation detail

For each individual, DNA was eluted from two extraction wells (40 µL per well) and
pooled, giving 80 µL of extract. Four separate 15-µL extract aliquots were used for
library preparation. Full-UDG libraries a, b, and c received 3, 2, and 1 µL of USER
enzyme, respectively. The non-UDG library used 3 µL of USER during loop-adaptor
cleavage but did not receive the dedicated full-UDG damage-removal treatment. All
four libraries underwent 12 indexing-PCR cycles.

## Read preprocessing, mapping, and duplicate removal

Adapters were removed and overlapping pairs were collapsed with AdapterRemoval v2.2.2. The operative settings were `--collapse`, `--trimns`, `--trimqualities`, `--minlength 30`, and `--minquality 25`. Three read sets were evaluated for each library:

1. collapsed reads;
2. collapsed plus `collapsed.truncated` reads; and
3. collapsed and `collapsed.truncated` reads plus unmerged paired-end reads.

Collapsed and `collapsed.truncated` reads were mapped as single-end data to hs37d5 with `bwa aln` followed by `bwa samse`. Unmerged `pair1.truncated` and `pair2.truncated` reads were mapped with `bwa aln` followed by `bwa sampe`. BAMs were converted, coordinate-sorted, merged across sequencing runs, and indexed with samtools.

For sets containing paired-end alignments, duplicate processing used name sorting, `samtools fixmate -m`, coordinate sorting, and `samtools markdup -r`. Single-end collapsed sets were coordinate-sorted and processed with `samtools markdup -s -r`. All three read sets were evaluated through PMD handling, but the reported imputation and IBD analyses used collapsed reads only. Adding unmerged reads contributed little coverage while reducing the endogenous-DNA proportion, which motivated this choice.

## PMD correction states

Each library was represented by an uncorrected BAM and correction-specific BAMs.
The token `raw` is retained only as an internal file/code identifier for these
standard-processed, uncorrected BAMs; it does not mean unprocessed sequencing data.

- **Terminal trimming:** trimBam soft-clipped 5 bases at each end of full-UDG reads and 5 or 10 bases at each end of non-UDG reads. Soft-clipped bases remain in the BAM but are not used for genotype-likelihood calculation.
- **mapDamage2 rescaling:** mapDamage2 was run with rescaling and merged-library modelling. The recorded flags were `--rescale` and `--merge-libraries`, with `--seq-length 5`, `--seq-length 10`, or the program default of 12 bases.
- **bamRefine masking:** bamRefine used the GRCh37-coordinate 1000 Genomes biallelic SNP list and PMD-length thresholds of 5 or 10 bases.

The default mapDamage2 state is called **Rescale-12** in the article and is stored
as `rescaled` in some archived paths and manifests. Not every generated state was
used for every endpoint. The exact endpoint-specific sets were:

| Endpoint | Full-UDG states | Non-UDG states |
|---|---|---|
| Alternate-read fractions and site retention | raw, Trim-5, Rescale-5, Rescale-12, bamRefine-5 | raw, Trim-5, Trim-10, Rescale-5, Rescale-10, Rescale-12, bamRefine-5, bamRefine-10 |
| Fixed-comparator imputation concordance | raw only, held fixed as comparator | all eight non-UDG states listed above |
| ancIBD and TKGWV2 state comparisons | raw, Trim-5, Rescale-5, bamRefine-5 | raw, Trim-5, Rescale-5, bamRefine-5 |
| KING-robust and IBS0 analysis of ORT15 and ORT16 | raw full-UDG pair only | not used |
| Matched-48 PCA and coverage-titration PCA | raw, Trim-5, Rescale-5, bamRefine-5 | all eight non-UDG states listed above |

mapDamage2 profiles were regenerated after correction as a diagnostic. mapDamage2
counts observed substitutions without weighting the revised base qualities, while
bamRefine changes only bases overlapping selected SNPs. Their post-correction curves
therefore need not flatten. Soft-clipping removes terminal bases from the plotted
alignment and visibly flattens the read-end curve.

For Supplementary Figure S4, the archived full-UDG profile tables came from merged
pre-deduplication collapsed BAMs, whereas the non-UDG profile tables came from
post-deduplication collapsed BAMs. The figure labels and caption state this distinction;
the diagnostic supports within-row comparisons among correction states rather than a
direct quantitative comparison between library types.

## Alternate-read fractions and site retention

Reference and alternate read counts were measured at biallelic 1000 Genomes sites with bcftools mpileup. Reads required MAPQ $\geq 30$ and base quality $\geq 20$; indels were excluded and the depth cap was 10,000. `FORMAT/DP` and `FORMAT/AD` were emitted, and REF/ALT alleles were aligned to the 1000 Genomes site panel before summarisation.

For a biallelic record, `FORMAT/AD` was interpreted as:

- `AD[0]`: reads supporting REF;
- `AD[1]`: reads supporting ALT; and
- analysed depth: `AD[0] + AD[1]`.

The site-level alternate-allele fraction was

\[
f_{ALT}=\frac{AD[1]}{AD[0]+AD[1]}.
\]

Two denominators were deliberately used:

- Transition and transversion means used only sites with `AD[0] + AD[1] >= 3`. Damage-associated transitions were C→T and G→A; transversions formed the comparison group.
- Site retention counted all covered panel sites (`AD[0] + AD[1] > 0`) before and after correction.

The retention denominator is therefore larger than the denominator for mean allele fractions. Reported fractions were calculated from `AD[0] + AD[1]`, rather than from a depth field that could include reads absent from the allele-depth counts.

For each correction state, the plotted transition–transversion difference is the mean
C→T/G→A alternate-read fraction minus the mean transversion fraction. Each correction state
retains its own callable site set; the analysis does not impose a fixed common-site
intersection. The difference is therefore a descriptive measure of preprocessing,
not a direct test of genuine variant loss.

## Imputation and fixed-comparator concordance

Autosomal genotype likelihoods were calculated at biallelic 1000 Genomes Phase 3 SNPs with bcftools 1.16/htslib 1.16, MAPQ $\geq 30$, and base quality $\geq 20$. Indels were excluded, no minimum read-depth filter was imposed, and the executed mpileup used the bcftools 1.16 default per-file maximum depth of 250. The 250-read setting is an upper cap, not a minimum. In particular, the combined REF+ALT depth $\geq 3$ rule used for the alternate-read-fraction endpoint was not applied to imputation; a site covered by one or two eligible reads could contribute genotype likelihoods. GLIMPSE v1.1.1 used the 2,504-sample Phase 3 SHAPEIT2 integrated v5b reference panel and GRCh37 genetic maps. Chromosomes were divided into chunks with a 2-Mb minimum window and 200-kb buffer, imputed with `GLIMPSE_phase`, ligated, and converted to phased GTs with `GLIMPSE_sample --solve`. The public wrapper states seed 15052011 explicitly; this is the fixed GLIMPSE v1.1.1 default used when the executed command omitted `--seed`. Because the sampled output contains GT but not GP, posterior genotype probabilities were copied from the ligated BCF, sample identifiers were normalised, and chromosomes 1–22 were concatenated. Phase accuracy and switch error were not evaluated separately.

The eight non-UDG states were imputed independently. For concordance, each was
compared with the same individual's uncorrected, independently imputed full-UDG
dataset (`raw` in code); no corrected full-UDG state was substituted between comparisons.
The full-UDG comparator is another low-coverage library from the same individual,
not an independent genotype truth. The concordance statistics therefore quantify
agreement between library preparations, and shared reference-panel errors can
increase apparent agreement without establishing accuracy.

Posterior dosage was reconstructed for every target VCF from genotype probabilities ordered as P(0/0), P(0/1), and P(1/1):

\[
DS=GP_1+2GP_2.
\]

The executed concordance inputs stored reconstructed dosage to four decimal places. Dosage \(r^2\) uses posterior dosage, whereas best-guess \(r^2\) uses the target maximum-posterior genotype derived from GP, not the sampled phased GT. The two statistics are therefore distinct.

GLIMPSE2_concordance v2.0.0 (commit 3bed6d9; release 2022-12-07) used the fixed full-UDG comparator, `--gt-val`, and the reference-frequency `AF` tag. The final run did not force target hard genotypes and did not bin by ALT frequency; neither `--gt-tar` nor `--use-alt-af` was used. No minimum depth was imposed on the fixed comparator, and no target-GP threshold was imposed. Frequency boundaries were 0, 0.001, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, and 0.50. Reported dosage and best-guess \(r^2\) values are count-weighted means across bins 1–8 (MAF $\geq 0.1\%$); bin 0 was excluded. NRD pooled all concordance-evaluable non-reference genotype comparisons across bins 0–8: all discordant genotypes divided by correctly matched heterozygous and alternate-homozygous genotypes plus all discordant genotypes. Correctly matched reference-homozygous pairs and sites without an eligible frequency bin do not enter this denominator.

## Paired leave-one-chromosome-out jackknife

The paired estimand was defined separately for each individual as

\[
\Delta NRD=NRD_{bamRefine\text{-}5}-NRD_{Rescale\text{-}5}.
\]

GLIMPSE2_concordance was rerun once per chromosome for the two correction methods in each individual, using the same comparator, filters, frequency panel, and options as the genome-wide run. Per-chromosome NRD numerators and denominators were extracted from the GLIMPSE confusion counts. Delete-one estimates \(d_i\) were formed by summing counts over the 21 retained autosomes, not by taking an unweighted mean of chromosome-specific percentages.

For \(n=22\) autosomes, the jackknife standard error was

\[
SE=\sqrt{\frac{n-1}{n}\sum_{i=1}^{n}(d_i-\bar d)^2}.
\]

The full-data difference was denoted \(\hat\theta\), and the bias-corrected
jackknife estimate was

\[
\hat\theta_{J}=n\hat\theta-(n-1)\bar d.
\]

The 95% interval was \(\hat\theta_J \pm t_{21}SE\), with \(t_{21}=2.080\).
The output reports both the direct full-data difference and the bias-corrected
jackknife estimate. Summed per-chromosome counts reproduced the genome-wide NRD
within 0.01 percentage points for all four individual-by-method datasets; the
largest observed deviation from the rounded table values was 0.004 percentage
points. The resulting intervals span zero and do not demonstrate statistical
equivalence.

## Final ancIBD v62 IBD1 analysis

The reported IBD analysis used ancIBD v0.7 on an exact allele-matched intersection of:

1. AADR v62 autosomal targets;
2. 1000 Genomes reference-frequency records; and
3. sites present across the GLIMPSE inputs.

The final intersection contained 1,100,313 SNPs. Genetic positions were supplied in Morgans from AADR v62. Native IBD1 calling began at 8 cM and retained ancIBD's internal post-processing. The final workflow was IBD1-only: IBD2 was not called, and no external gap merging was applied.

Summaries applied strict length cutoffs of >8, >12, >16, and >20 cM. The primary filter retained segments strictly longer than 12 cM with marker density strictly greater than 220 SNP/cM. Both total retained IBD1 length and segment count were recorded. The IBD1 proportion used the summed first-to-last callable marker span across chromosomes 1–22, 3,538.8511 cM, as its denominator.

For each individual, the ancIBD state axis contained eight inputs in a fixed order:
four full-UDG states (raw, Trim-5, Rescale-5, and bamRefine-5) followed by the same
four non-UDG states. The Cartesian product of the eight ORT15 states and eight
ORT16 states produced the complete ordered 8×8 matrix of 64 cells. It includes
same-state, cross-state, same-library-type, and mixed-library-type comparisons;
reflected mixed-library cells retain different individual-by-library assignments
and are not technical replicates.

Two prespecified subsets were used for the main summaries:

- The asymmetric subset contains eight comparisons: one individual's full-UDG
  state remained raw while the other individual's non-UDG state varied across
  raw, Trim-5, Rescale-5, and bamRefine-5, in both directions.
- The same-method subset contains 16 comparisons: each of the four state labels
  was applied to full-UDG×full-UDG, non-UDG×non-UDG, and the two directional
  mixed-library pairings.

The two raw mixed-library cells occur in both subsets; their union therefore
contains 22 distinct cells. The other 42 cells are cross-state comparisons shown
for matrix context rather than as one-factor or same-method contrasts. Larger
inferred IBD1 length measures greater inferred sharing for that comparison; it is
not an accuracy ranking in the absence of independent truth and unrelated controls.

## READv2 pseudo-haploid calling and normalisation

The primary READv2 result used an archived cemetery-wide pseudo-haploid PLINK
callset at 1240k targets and was not a correction-state sweep. The genotype input
was generated from ORT1–ORT27 and included chromosomes 1–22, X, and Y. READv2
returned 210 pairwise rows among 21 retained individuals. For ORT15 and ORT16,
the executed upstream merge included only full-UDG libraries a, b, and c; the
merged BAM was then soft-clipped by 6 bp at each read end. The run used `samtools
mpileup -B -q30 -Q30 -R` and pileupCaller `--randomHaploid` to sample one observed
allele per covered target. Normalisation used the median pairwise \(P_0\) as the
unrelated baseline.

The cohort was predominantly, but not entirely, unrelated; it included ORT15,
ORT16, and other inferred relatives. No leave-one-out or relative-excluded
sensitivity analysis was run, so the result should not be described as insensitive
to cohort composition. The primary ORT15–ORT16 result used 169,172 overlapping
targets (165,350 autosomal, 3,684 X-chromosome, and 138 Y-chromosome targets) and
returned KC 0.230967 with the parent–offspring subtype. A later autosome-only run
using ORT6–ORT27 returned KC 0.231592 at 157,451 overlapping SNPs and the same
first-degree parent–offspring classification. Only these two ascertained-panel
READv2 analyses are included as reported evidence.

## KING-robust kinship and IBS0

KING-robust kinship and IBS0 were calculated from the same merged, GLIMPSE-imputed
diploid genotypes. The inputs were the uncorrected full-UDG collapsed-only
VCFs `ORT15_fu.1.imputed.allchr.vcf.gz` and
`ORT16_fu.1.imputed.allchr.vcf.gz`; `.1` denotes the collapsed-only read set.
No PMD-corrected VCF was substituted for this calculation. IBS0 sites are
opposing homozygotes: one individual is 0/0 and the other 1/1. Genome-wide analysis
was restricted to 5,272,558 variants with reference-panel allele frequency \(p\)
between 0.1 and 0.9, with \(q=1-p\).

The site-specific unrelated opposing-homozygote probability is \(2p^2q^2\). The sample-specific expectations were therefore

\[
E_{unrelated}=\operatorname{mean}(2p^2q^2)
\]

and

\[
E_{full\ sibling}=\frac{1}{4}E_{unrelated},
\]

because full siblings are IBD0 across approximately one quarter of the genome. Under a parent–offspring model, the theoretical IBS0 rate is zero except for genotype error. KING-robust kinship was computed as

\[
\frac{N_{het,het}-2N_{IBS0}}{N_{het,1}+N_{het,2}}.
\]

To evaluate residual imputation error, the IBS0 analysis was repeated on chromosome 1 after retaining sites with posterior genotype probability `GP >= 0.99` in both individuals. This left 324,285 sites; the expected rates were recomputed from the allele frequencies of that same subset rather than copied from the genome-wide set.

A same-individual read-set control compared ORT15 full-UDG collapsed-only data with the nested collapsed-plus-truncated read set on chromosome 1 at reference-panel allele frequency 0.1–0.9. It yielded 5 IBS0 calls among 408,169 filtered sites, or $1.23\times10^{-5}$ per all filtered site. The alternative value 5/267,326 (0.00001870, rounded to 0.00002) uses only jointly homozygous sites as the denominator and is not the per-all-site rate used for ORT15 and ORT16.

## Standalone mitochondrial contamination analysis

Mitochondrial contamination was assessed on collapsed, deduplicated non-UDG BAMs with the mtCont component of schmutzi v1.5.7. These libraries show asymmetric damage: reduced 5$'$ C→T and strong 3$'$ G→A. Accordingly, empirical double-stranded 5$'$ and 3$'$ damage profiles were supplied separately instead of applying a symmetric profile. endoCaller supplied endogenous-base likelihoods and was run under a single-contaminant-haplotype assumption. mtCont compared the data with the schmutzi Eurasian mitochondrial frequency panel.

This was a standalone, one-pass mtCont analysis, not the complete iterative schmutzi workflow. The reported values are the program's mtCont estimates and intervals. An unmodelled-damage mtCont run and A/T-site minor-allele summaries were sensitivity/QC analyses and are not the reported contamination estimator.

## Mitochondrial variant calling and haplogroup assignment

Mitochondrial variants were called from each full-UDG `MT` BAM with bcftools v1.16/htslib v1.16 as follows:

```text
bcftools mpileup -r MT -f hs37d5.fa -B -q30 -Q30 -Ou FULL_UDG_MTY.bam |
bcftools call --ploidy 1 -m
```

Thus, the executed call set used MAPQ $\geq 30$ and base quality $\geq 30$ with BAQ disabled (`-B`) and haploid calling. Direct inspection of the complete VCFs confirms that no global depth $\geq 3$ filter was applied: 52 ORT15 records and 173 ORT16 records have `INFO/DP` of 1 or 2.

The low-depth records did not enter the called SNP set used for the haplogroup interpretation. Each individual has 37 called single-nucleotide ALT variants, and the minimum `INFO/DP` among those variants is 5. Consequently, applying a post hoc depth $\geq 3$ criterion to the 37 called SNPs would remove none of them and would not change that haplogroup variant set. This observation does not establish that a depth-filtered whole-mitogenome consensus was executed; it establishes only that all called haplogroup SNPs already exceed the stated depth threshold.

HaploGrep 3 v3.3.2 was subsequently used with PhyloTree 17 Forensic Update (`phylotree-fu-rcrs@1.3`), producing D4o1 assignments with reported quality scores 0.891 and 0.911. For the between-individual comparison, called single-nucleotide ALT records were keyed by position, REF, and ALT. Each individual had 37 such SNPs; 36 were shared at the same position and allele. The individual-specific SNPs were ORT15 A4215G and ORT16 C7028T, both at depth 5. The additional ORT15 record at position 3106 (`CN` to `C`) was treated separately as an rCRS reference-ambiguity record rather than as one of the 37 SNPs.

## PCA projection and marker trace

Twenty-four correction-state datasets were analysed: 12 per individual. The four
full-UDG states were raw, Trim-5, Rescale-5, and bamRefine-5. The eight non-UDG
states were raw, Trim-5, Trim-10, Rescale-5, Rescale-10, Rescale-12,
bamRefine-5, and bamRefine-10. Each was represented both as GLIMPSE-imputed
genotypes and as non-imputed pseudo-haploid calls produced by sampling one observed
allele at each covered 1240k site with pileupCaller `--randomHaploid`. This yielded
48 projected query points on a shared reference basis.

Reference axes were constructed in smartpca (EIGENSOFT build 18140) from AADR v62.0. The reference included 985 present-day HGDP and SGDP individuals from 154 populations; ancient and archaic individuals were excluded. After allele alignment, 1,134,702 autosomal markers were shared between reference and query panels; 1,113,348 entered smartpca after internal filtering.

ORT samples were projected without influencing the axes using:

```text
lsqproject: YES
numoutlieriter: 0
```

No LD pruning was applied. PC1–PC4 explained 6.45%, 4.37%, 1.53%, and 1.04% of total variance. Each percentage was calculated as the component eigenvalue divided by the full trace of all 984 non-zero components, not by the sum of the displayed components. Correction-associated spread was the maximum pairwise Euclidean PC distance within each individual-by-library-by-representation group and was compared descriptively with the scatter among Yakut reference individuals.

## Coverage-titration PCA

The coverage titration was a separate, non-imputed pseudo-haploid smartpca run. It
used the same four full-UDG and eight non-UDG correction states per individual as
the matched-48 PCA. Current autosomal mean depth was measured from chromosomes
1–22 with `samtools coverage`. Each BAM was first restricted to the 1240k target
intervals; uniform subsampling of that restricted BAM is equivalent for calls at
those targets and avoided rescanning the full BAM for every replicate.

The discrete target depths were 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, and 0.20×.
For each target below the source depth, the retained-read fraction was target depth
divided by measured source depth. Ten BAM subsamples were generated with
`samtools view -s`, using integer downsampling seeds 1–10 and the six-digit retained
fraction in samtools' `seed.fraction` argument. One undownsampled full-coverage
anchor was also included per correction state. ORT15 non-UDG Trim-10 had a measured
source depth of 0.172982× and therefore could not produce the ten 0.20× replicates;
this was the only unattainable state-by-target cell. The final projection contained
1,694 ORT rows: 842 for ORT15 and 852 for ORT16.

For each coverage level, `samtools mpileup -B -q30 -Q30 -R` was run at the 1240k
targets, followed by pileupCaller `--randomHaploid`; no imputation was performed.
The archived pileupCaller command does not pass an explicit random seed, so the
recorded replicate seeds refer to BAM subsampling. All coverage levels were combined
as PLINK data and projected in one smartpca run onto 985 present-day HGDP/SGDP AADR
v62 individuals from 154 populations. Marker IDs were harmonised as chromosome and
position, strand-conflict sites were excluded on merge, and ORT rows were labelled
as projected. The separate run used `lsqproject: YES`, `numoutlieriter: 0`, and
`numoutevec: 10`, with no LD-pruning step. Its coordinates are internally comparable
within the titration but are not the same coordinate set as the matched-48 PCA.
