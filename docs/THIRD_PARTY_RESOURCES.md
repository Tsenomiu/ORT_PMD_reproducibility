# Third-party resources

Reference panels, reference genomes, and third-party software are not redistributed
in this repository. Obtain them from their original providers and follow the
applicable terms of use.

## Reference data

| Resource | Version/build used | Role |
|---|---|---|
| hs37d5 | GRCh37/hs37d5 | Read mapping and pileup reference |
| 1000 Genomes Project | Phase 3, GRCh37 | SNP sites, GLIMPSE reference haplotypes, allele frequencies, genetic maps, and bamRefine site list |
| Allen Ancient DNA Resource | AADR v62.0, 1240k | ancIBD marker positions, PCA reference data, and pseudo-haploid target sites |
| HGDP and SGDP through AADR | 985 individuals from 154 populations | PCA reference axes |
| HapMapChrX | Version distributed with the ANGSD analysis environment | Male-X contamination estimation |
| schmutzi Eurasian panel | schmutzi v1.5.7 distribution | mtCont frequency comparison |
| Revised Cambridge Reference Sequence | NC_012920 | Mitochondrial variant calling |
| PhyloTree 17 Forensic Update | `phylotree-fu-rcrs@1.3` | HaploGrep 3 haplogroup assignment |

Record the provider version, local filename, and SHA-256 checksum for every resource
used in a rerun. Do not commit completed local configurations or third-party data.

## Software

The workflows use AdapterRemoval, BWA, samtools/bcftools, trimBam, mapDamage2,
bamRefine, GLIMPSE, GLIMPSE2_concordance, ancIBD, READv2, pileupCaller, PLINK,
ANGSD, schmutzi, HaploGrep 3, and smartpca/EIGENSOFT. Exact study versions are listed
in [`../software_versions.tsv`](../software_versions.tsv).

Install each program from its official project and cite it as requested by its
authors. The repository does not vendor third-party executables or source trees.

TKGWV2 contributed a result reported in the associated article, but a public runner
is not included because the locally patched conversion step is not available as a
verified portable script.

## Mitochondrial phylogeny citations

Use of `phylotree-fu-rcrs@1.3` requests citation of both:

- Weissensteiner H, et al. *HaploGrep 2: mitochondrial haplogroup classification in
  the era of high-throughput sequencing.* Nucleic Acids Research.
  2016;44(W1):W58–W63. doi:10.1093/nar/gkw233.
- Dür A, Huber N, Parson W. *Fine-Tuning Phylogenetic Alignment and Haplogrouping of
  mtDNA Sequences.* International Journal of Molecular Sciences.
  2021;22(11):5747. doi:10.3390/ijms22115747.

These citations do not grant permission to redistribute the phylogeny.
