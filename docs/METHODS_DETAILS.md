# Additional methods

This document provides computational details that complement the associated
article. Executable commands and input schemas are organised by analysis stage in
[`../workflows/`](../workflows/).

## Read preprocessing, mapping, and duplicate removal

Adapters were removed and overlapping pairs were collapsed with AdapterRemoval v2.2.2. The operative settings were `--collapse`, `--trimns`, `--trimqualities`, `--minlength 30`, and `--minquality 25`. Three read sets were evaluated for each library:

1. collapsed reads;
2. collapsed plus `collapsed.truncated` reads; and
3. collapsed and `collapsed.truncated` reads plus unmerged paired-end reads.

Collapsed and `collapsed.truncated` reads were mapped as single-end data to hs37d5 with `bwa aln` followed by `bwa samse`. Unmerged `pair1.truncated` and `pair2.truncated` reads were mapped with `bwa aln` followed by `bwa sampe`. BAMs were converted, coordinate-sorted, merged across sequencing runs, and indexed with samtools.

For sets containing paired-end alignments, duplicate processing used name sorting, `samtools fixmate -m`, coordinate sorting, and `samtools markdup -r`. Single-end collapsed sets were coordinate-sorted and processed with `samtools markdup -s -r`. All three read sets were evaluated through PMD handling, but the reported imputation and IBD analyses used collapsed reads only. Adding unmerged reads contributed little coverage while reducing the endogenous-DNA proportion, which motivated this choice.

## PMD-handling treatments

Each library was represented by an uncorrected BAM and treatment-specific BAMs.

- **Terminal trimming:** trimBam soft-clipped 5 bases at each end of full-UDG reads and 5 or 10 bases at each end of non-UDG reads. Soft-clipped bases remain in the BAM but are not used for genotype-likelihood calculation.
- **mapDamage2 rescaling:** mapDamage2 was run with rescaling and merged-library modelling. The recorded flags were `--rescale` and `--merge-libraries`, with `--seq-length 5`, `--seq-length 10`, or the program default of 12 bases.
- **bamRefine masking:** bamRefine used the GRCh37-coordinate 1000 Genomes biallelic SNP list and PMD-length thresholds of 5 or 10 bases.

mapDamage2 profiles were regenerated after treatment as a diagnostic. Rescaling and bamRefine masking alter the contribution of bases to downstream likelihoods; they do not remove the physical terminal mismatches, so their post-treatment misincorporation curves need not flatten in the same way as trimming curves.

## Reference and alternative allele counts

Reference and alternative read counts were measured at biallelic 1000 Genomes sites with bcftools mpileup. Reads required MAPQ $\geq 30$ and base quality $\geq 20$; indels were excluded and the depth cap was 10,000. `FORMAT/DP` and `FORMAT/AD` were emitted, and REF/ALT alleles were aligned to the 1000 Genomes site panel before summarisation.

For a biallelic record, `FORMAT/AD` was interpreted as:

- `AD[0]`: reads supporting REF;
- `AD[1]`: reads supporting ALT; and
- analysed depth: `AD[0] + AD[1]`.

The site-level alternative-allele fraction was

\[
f_{ALT}=\frac{AD[1]}{AD[0]+AD[1]}.
\]

Two denominators were deliberately used:

- Transition and transversion means used only sites with `AD[0] + AD[1] >= 3`. Damage-associated transitions were C→T and G→A; transversions formed the comparison group.
- Site retention counted all covered panel sites (`AD[0] + AD[1] > 0`) before and after treatment.

The retention denominator is therefore larger than the denominator for mean allele fractions. Reported fractions were calculated from `AD[0] + AD[1]`, rather than from a depth field that could include reads absent from the allele-depth counts.

## Imputation and fixed-comparator concordance

Genotype likelihoods were calculated on chromosomes 1–22 at 1000 Genomes Phase 3 sites with MAPQ $\geq 30$ and base quality $\geq 20$. GLIMPSE v1 used the Phase 3 reference panel and genetic maps for chunked imputation and phasing, followed by chromosome-level concatenation.

For concordance, each individual's unmodified, imputed full-UDG dataset was held fixed as the common matched comparator for every non-UDG treatment. It is not an independent genotype truth.

Posterior dosage was reconstructed for every target VCF from genotype probabilities ordered as P(0/0), P(0/1), and P(1/1):

\[
DS=GP_1+2GP_2.
\]

Dosage \(r^2\) uses posterior-mean dosage, whereas best-guess \(r^2\) uses the hard genotype call. Reconstructing `DS` from `GP` therefore keeps the two statistics distinct.

GLIMPSE2_concordance used the fixed full-UDG comparator, `--gt-val`, and the reference-frequency `AF` tag. The final run did not force target hard genotypes and did not bin by ALT frequency; neither `--gt-tar` nor `--use-alt-af` was used. Frequency boundaries were 0, 0.001, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, and 0.50. Reported dosage and best-guess \(r^2\) values are count-weighted means across bins 1–8 (MAF $\geq 0.1\%$); bin 0 was excluded. NRD was calculated across all evaluated sites.

## Paired leave-one-chromosome-out jackknife

The paired estimand was defined separately for each individual as

\[
\Delta NRD=NRD_{bamRefine\text{-}5}-NRD_{Rescale\text{-}5}.
\]

GLIMPSE2_concordance was rerun once per chromosome for the two treatments in each individual, using the same comparator, filters, frequency panel, and options as the genome-wide run. Per-chromosome NRD numerators and denominators were extracted from the GLIMPSE confusion counts. Delete-one estimates \(d_i\) were formed by summing counts over the 21 retained autosomes, not by taking an unweighted mean of chromosome-specific percentages.

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
within 0.01 percentage points for all four individual-by-treatment datasets; the
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

Two comparison designs were retained:

- In the asymmetric design, one individual's full-UDG library remained raw while only the other individual's non-UDG treatment varied, with both directions evaluated.
- In the same-treatment design, the same treatment was applied to both members of the full-UDG×full-UDG, non-UDG×non-UDG, and two directional mixed-library pairings.

The complete ordered matrix contains 64 comparisons. Directional mixed pairs were not collapsed because the individual, library type, and coverage differ by side. Larger inferred IBD1 length measures greater inferred sharing for that comparison; it is not an accuracy ranking in the absence of independent truth and unrelated controls.

## READv2 pseudo-haploid calling and normalisation

READv2 used pseudo-haploid PLINK genotypes at 1240k sites. The 22-individual Oortsog Ovoo comparison cohort comprised ORT6–ORT27. Reads were piled up at the 1240k panel and pileupCaller sampled one observed allele per covered site with `--randomHaploid`. Autosomes were retained for READv2. Normalisation used the cohort median pairwise \(P_0\) as the unrelated baseline.

The cohort was predominantly, but not entirely, unrelated; it included ORT15, ORT16, and other inferred relatives. No leave-one-out or relative-excluded sensitivity analysis was run, so the result should not be described as insensitive to cohort composition. The reported ORT15–ORT16 1240k result used 169,172 overlapping SNPs. The first-degree subtype was stable on this ascertained panel; all-1000-Genomes-site analyses were treated as sensitivity analyses because subtype calls changed under rescaling.

## KING-robust kinship and IBS0

KING-robust kinship and IBS0 were calculated from the same merged, GLIMPSE-imputed diploid genotypes. IBS0 sites are opposing homozygotes: one individual is 0/0 and the other 1/1. Genome-wide analysis was restricted to 5,272,558 variants with reference-panel allele frequency \(p\) between 0.1 and 0.9, with \(q=1-p\).

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

A same-individual read-set control compared ORT15 full-UDG collapsed-only data with the nested collapsed-plus-truncated read set on chromosome 1 at reference-panel allele frequency 0.1–0.9. It yielded 5 IBS0 calls among 408,169 filtered sites, or $1.23\times10^{-5}$ per all filtered site. The alternative value 5/267,326 (0.00001870, rounded to 0.00002) uses only jointly homozygous sites as the denominator and is not the per-all-site rate used for the focal pair.

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

HaploGrep 3 v3.3.2 was subsequently used with PhyloTree 17 Forensic Update (`phylotree-fu-rcrs@1.3`), producing D4o1 assignments with reported quality scores 0.891 and 0.911.

## PCA projection and marker trace

Twenty-four treatment datasets were analysed: 12 per individual, comprising four full-UDG and eight non-UDG treatments. Each was represented both as GLIMPSE-imputed genotypes and as non-imputed pseudo-haploid calls produced by sampling one observed allele at each covered 1240k site with pileupCaller `--randomHaploid`. This yielded 48 projected query points on a shared reference basis.

Reference axes were constructed in smartpca (EIGENSOFT build 18140) from AADR v62.0. The reference included 985 present-day HGDP and SGDP individuals from 154 populations; ancient and archaic individuals were excluded. After allele alignment, 1,134,702 autosomal markers were shared between reference and query panels; 1,113,348 entered smartpca after internal filtering.

ORT samples were projected without influencing the axes using:

```text
lsqproject: YES
numoutlieriter: 0
```

No LD pruning was applied. PC1–PC4 explained 6.45%, 4.37%, 1.53%, and 1.04% of total variance. Each percentage was calculated as the component eigenvalue divided by the full trace of all 984 non-zero components, not by the sum of the displayed components. Treatment-associated spread was the maximum pairwise Euclidean PC distance within each individual-by-library-by-representation group and was compared descriptively with the scatter among Yakut reference individuals.
