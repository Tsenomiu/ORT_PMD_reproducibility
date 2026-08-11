# Third-party resources

Reference panels, reference genomes, and third-party software are not redistributed
in this repository. Obtain them from their original providers and follow the
applicable terms of use.

## Reference data

| Resource | Version/build used | Role |
|---|---|---|
| hs37d5 | GRCh37/hs37d5 | Read mapping and pileup reference |
| 1000 Genomes Project | Phase 3, GRCh37; SHAPEIT2 integrated v5b for GLIMPSE | SNP sites, GLIMPSE reference haplotypes, allele frequencies, genetic maps, and bamRefine site list |
| Allen Ancient DNA Resource | AADR v62.0, 1240k | ancIBD marker positions, PCA reference data, and pseudo-haploid target sites |
| HGDP and SGDP through AADR | 985 individuals from 154 populations | PCA reference axes |
| HapMapChrX | Exact local build/checksum not recorded; version distributed with the ANGSD analysis environment | Male-X contamination estimation |
| schmutzi Eurasian panel | schmutzi v1.5.7 distribution | mtCont frequency comparison |
| Revised Cambridge Reference Sequence | NC_012920 | Mitochondrial variant calling |
| PhyloTree 17 Forensic Update | `phylotree-fu-rcrs@1.3` | HaploGrep 3 haplogroup assignment |

Record the provider version, local filename, and SHA-256 checksum for every resource
used in a rerun. Do not commit completed local configurations or third-party data.

## Software

The workflows use AdapterRemoval, BWA, samtools/bcftools, trimBam, mapDamage2,
bamRefine, GLIMPSE, GLIMPSE2_concordance, ancIBD, READv2, pileupCaller, PLINK 1.9,
PLINK 2,
ANGSD, schmutzi, HaploGrep 3, Qualimap, and smartpca/EIGENSOFT. Recorded study
versions, including explicit unknowns, are listed in
[`../software_versions.tsv`](../software_versions.tsv).

Install each program from its official project and cite it as requested by its
authors. The repository does not vendor third-party executables or source trees.

Primary citations added to the associated article for version-specific tools
include Jun et al. (2015; doi:10.1101/gr.176552.114) for BamUtil/trimBam,
Rubinacci et al. (2023; doi:10.1038/s41588-023-01438-3) for GLIMPSE2,
Patterson et al. (2006; doi:10.1371/journal.pgen.0020190) for the smartpca
eigenanalysis method, and Schiffels (2026; doi:10.21105/joss.08634) for
pileupCaller.

TKGWV2 v1.0b is pinned at upstream commit
`c8638d47d3143b82ec969259df66e16f63b2b0ae`. The release candidate includes the
exact recovered four-line compatibility patch and a processed-input runner
reconstructed from successful logs. TKGWV2 and the patch are covered by
GPL-2.0-only; see `../workflows/08_kinship/tkgwv2/THIRD_PARTY_NOTICE.md`.

## Mitochondrial reference, software and phylogeny citations

The mitochondrial analysis uses rCRS, HaploGrep 3 and
`phylotree-fu-rcrs@1.3`. Cite:

- Andrews RM, et al. *Reanalysis and revision of the Cambridge reference
  sequence for human mitochondrial DNA.* Nature Genetics. 1999;23:147.
  doi:10.1038/13779.
- Schönherr S, Weissensteiner H, Kronenberg F, Forer L. *Haplogrep 3—an
  interactive haplogroup classification and analysis platform.* Nucleic Acids
  Research. 2023;51(W1):W263–W268. doi:10.1093/nar/gkad284.

The selected phylogeny additionally requests citation of both:

- Weissensteiner H, et al. *HaploGrep 2: mitochondrial haplogroup classification in
  the era of high-throughput sequencing.* Nucleic Acids Research.
  2016;44(W1):W58–W63. doi:10.1093/nar/gkw233.
- Dür A, Huber N, Parson W. *Fine-Tuning Phylogenetic Alignment and Haplogrouping of
  mtDNA Sequences.* International Journal of Molecular Sciences.
  2021;22(11):5747. doi:10.3390/ijms22115747.

These citations do not grant permission to redistribute the phylogeny.
