# Reproducibility limitations

This repository distinguishes verification of the reported values from a
complete reconstruction starting with raw sequencing reads.

## READv2 boundary

The exact original primary BAM-to-genotype command and random state used for
pseudo-haploid sampling were not preserved. The derived genotype inputs,
downstream READv2 commands, normalization procedure and reported outputs are
provided to permit verification of the published calculations. The primary and
autosome-only 1240k analyses reproduce the reported parent–offspring result; the
transversion-only analysis supports first-degree relatedness but does not assign
a parent–offspring subtype.

The compact public package verifies the reported normalization arithmetic and
classifications. It does not claim byte-identical regeneration of the primary
pseudo-haploid genotype calls from BAM files.

## Raw and controlled data

Raw FASTQ, BAM, large VCF, PLINK genotype matrices, reference genomes, and
licensed reference panels are excluded from GitHub. Raw study reads are
registered under ENA study accession PRJEB112497. The ENA study remains private;
the accession should be treated as a public raw-data route only after ENA release
has been verified.

Controlled archaeological coordinates and raw micro-CT volumes are also not
distributed in this repository. Supplementary Figure S2 and Supplementary
Figure S3 therefore retain documented external-input boundaries.

## TKGWV2 boundary

All 16 reported TKGWV2 comparisons were reproduced exactly from the retained
processed PED/MAP boundary, including HRC 0.2356 and 4,517,214 used SNPs for the
full-UDG × full-UDG bamRefine-5 comparison. The BAM-to-PED pseudo-haploid step
used unseeded sampling, so deterministic reproduction is claimed only from the
processed-input boundary. Large PED/MAP files and the 1000 Genomes EAS frequency
resource are identified by checksum but are not distributed here.

## PCA boundary

The recovered workflow covers TEST query preparation, the authoritative RC
smartpca projection, method-distance calculations, and plotting. Public
verification uses ORT query-only coordinates and population-level reference
aggregates. Individual-level AADR/HGDP/SGDP reference coordinates and licensed
reference genotypes are excluded.

The archived 150-dpi PCA raster checks are byte-exact in the canonical Poppler
25.06.0 / TeX Live 2024 environment. The failing v1.1.0 GitHub Actions runs used
Poppler 24.02.0, which changed edge anti-aliasing and PNG bytes while leaving the
numerical outputs and plotted content unchanged. Local validation therefore
remains `canonical-strict`; CI explicitly uses `ci-portable`, which retains all 17
numerical and structural checks and requires the expected RGB PNG format,
dimensions, and documented visual-equivalence thresholds. Missing or visibly
different figures fail both modes.

## Archaeological report

The unpublished excavation report is not a public dataset and is excluded from
this repository. The repository does not reproduce its text or site-plan
source material. Figure 1 is an author-created workflow schematic and does not
require excavation-report adaptation credit.
