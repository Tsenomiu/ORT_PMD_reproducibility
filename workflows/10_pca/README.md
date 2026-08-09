# 10 Population genetic context

The PCA contains 24 GLIMPSE-imputed queries, the same 24 datasets represented as
random pseudo-haploid 1240k calls, and 985 present-day HGDP/SGDP individuals from
AADR v62. The two query representations are projected together onto shared reference
axes and do not influence those axes.

`run_matched48_projection.sh` uses PLINK 2 to restrict the inputs to autosomes and
normalise variant identifiers. It then uses PLINK 1.9 for marker extraction and
dataset merging, followed by EIGENSTRAT conversion, population labelling, and
smartpca projection of the three input PLINK datasets.

Required executables are PLINK 1.9, PLINK 2, EIGENSOFT `convertf`, EIGENSOFT
`smartpca`, and Python 3.

Final settings:

```text
smartpca / EIGENSOFT 18140
lsqproject: YES
numoutlieriter: 0
altnormstyle: NO
numoutevec: 10
```

No LD pruning is applied. The aligned intersection contains 1,134,702 sites;
smartpca uses 1,113,348 after internal filtering. Variance percentages use the full
trace of 984 non-zero components.

The compact coordinate check can be run independently:

```bash
python3 workflows/10_pca/verify_pca.py \
  data/summary/pca/pca_ort_queries_matched48.evec \
  data/summary/pca/pca_aadr_matched48.eval \
  data/summary/pca/query_manifest_matched48.tsv \
  data/summary/pca/reference_population_aggregates.tsv \
  work/pca_recomputed_metrics.csv
cmp data/summary/pca/pca_metrics_matched48.csv \
  work/pca_recomputed_metrics.csv
```

See [`QUERY_INPUTS.md`](QUERY_INPUTS.md) for the required query-data structure.
