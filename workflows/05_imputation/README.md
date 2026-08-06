# 05 Imputation and phasing

`run_glimpse1_chromosome.sh` performs the chromosome-level GLIMPSE v1 workflow:

1. calculate genotype likelihoods at Phase 3 sites with MAPQ 30 and base quality 20;
2. create chunks with a 2,000,000-bp window and 200,000-bp buffer;
3. phase against the Phase 3 reference panel and genetic map;
4. ligate chromosome chunks; and
5. sample phased genotypes with `GLIMPSE_sample --solve`.

Run the script for each treatment and chromosome 1–22. Input BAMs and 1000
Genomes/GLIMPSE resources are supplied through `workflows/config.sh`.

```bash
for chr in $(seq 1 22); do
  workflows/05_imputation/run_glimpse1_chromosome.sh \
    ORT15_nu_raw /path/to/input.bam "$chr" work/glimpse
done
```
