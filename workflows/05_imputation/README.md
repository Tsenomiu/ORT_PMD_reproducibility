# 05 Imputation and phased-genotype generation

`run_glimpse1_chromosome.sh` performs the chromosome-level GLIMPSE v1 workflow:

1. calculate genotype likelihoods at biallelic Phase 3 SNPs with MAPQ 30,
   base quality 20, no minimum depth, and the executed bcftools 1.16 maximum
   per-file depth of 250;
2. create chunks with a 2,000,000-bp window and 200,000-bp buffer;
3. phase against the Phase 3 reference panel and genetic map with seed
   15052011 (the fixed GLIMPSE v1.1.1 default, made explicit here);
4. ligate chromosome chunks; and
5. select the most likely phased haplotype pair with `GLIMPSE_sample --solve`;
6. restore FORMAT/GP from the ligated BCF and normalise the sample identifier; and
7. validate phased GT and GP fields and concatenate the 22 autosomes with
   `concat_glimpse1_autosomes.sh`.

The final GT field is phased, but this study did not evaluate switch error or phase
accuracy. FORMAT/GP is retained because the fixed-comparator concordance stage uses
target posterior probabilities and reconstructs dosage from them.

Run the script for each treatment and chromosome 1–22. Input BAMs and 1000
Genomes/GLIMPSE resources are supplied through `workflows/config.sh`.

```bash
for chr in $(seq 1 22); do
  workflows/05_imputation/run_glimpse1_chromosome.sh \
    ORT15_nu_raw ORT15 /path/to/input.bam "$chr" work/glimpse
done

workflows/05_imputation/concat_glimpse1_autosomes.sh \
  ORT15_nu_raw work/glimpse
```
