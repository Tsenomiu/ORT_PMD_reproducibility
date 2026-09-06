# Analysis workflows

The workflows are ordered by analysis stage. Each directory documents required
software, inputs, outputs, and any available executable example.

1. `01_sequence_processing`: adapter removal, mapping, merging, and duplicate removal
2. `02_contamination`: nuclear and mitochondrial contamination estimates
3. `03_pmd_correction`: trimBam, mapDamage2, and bamRefine treatments
4. `04_alt_fraction`: reference/alternate read counts, site retention, and DeltaALT uncertainty
5. `05_imputation`: GLIMPSE genotype likelihoods, imputation, and phased-genotype generation
6. `06_concordance`: fixed-comparator concordance, chromosome jackknife, and Table 1 generation
7. `07_ancibd`: AADR-v62 ancIBD analysis
8. `08_kinship`: TKGWV2, READv2, KING, and IBS0 analyses
9. `09_mtdna`: corrected collapsed-only mitochondrial calling, call-set comparison, and haplogroup assignment
10. `10_pca`: matched imputed and pseudo-haploid PCA preparation, projection, and recovery validation

## Configuration

```bash
cp workflows/config.example.sh workflows/config.sh
# Edit config.sh for the local installation and data locations.
export ORT_CONFIG="$PWD/workflows/config.sh"
source "$ORT_CONFIG"
```

The example configuration and manifests contain placeholders only. Raw FASTQ, BAM,
VCF, reference-panel, and genotype files are not stored in this repository.

## Scope

Many stages require server-scale compute and externally obtained reference resources.
The scripts are portable entry points, not a bundled computing environment. Follow
the README in each directory and obtain the resources listed in
[`../docs/THIRD_PARTY_RESOURCES.md`](../docs/THIRD_PARTY_RESOURCES.md).

Coordinates use GRCh37/hs37d5 except for mitochondrial rCRS NC_012920. Main
imputation, IBD and mitochondrial analyses use collapsed-only reads. Alternative
collapsed-plus-truncated and all-read states are nested supersets and must not be
merged. The full-UDG imputed dataset is a matched comparator, not an independent
genotype truth.
