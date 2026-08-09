# 09 Mitochondrial haplogroup assignment

`run_haplogroup.sh` merges the three deduplicated full-UDG sublibrary BAMs for an
individual, calls haploid mitochondrial variants against hs37d5/rCRS, reheaders the
sample, and runs HaploGrep 3.

The bcftools v1.16 call uses:

```text
bcftools mpileup -r MT -f hs37d5.fa -B -q30 -Q30 -Ou FULL_UDG.bam |
bcftools call --ploidy 1 -m
```

Haplogroup classification uses HaploGrep 3 v3.3.2 with
`phylotree-fu-rcrs@1.3`. Cite Andrews et al. (1999) for rCRS and Schönherr et
al. (2023) for HaploGrep 3. Also cite Weissensteiner et al. (2016) and Dür et
al. (2021), as requested by the phylogeny provider.

Input BAMs, called VCFs, and mitochondrial consensus data are not distributed.
