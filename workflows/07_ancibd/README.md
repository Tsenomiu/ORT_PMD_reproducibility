# 07 ancIBD v62 analysis

The workflow uses ancIBD v0.7 and an exact allele-matched intersection of AADR
v62, 1000 Genomes allele frequencies, and the 16 treatment VCFs.

Key settings:

- expected autosomal intersection: 1,100,313 SNPs;
- genetic coordinates supplied in Morgans;
- IBD1 only, with native calling beginning at 8 cM;
- ancIBD native post-processing and no external gap merge;
- summary thresholds strictly greater than 8, 12, 16, and 20 cM; and
- primary filter strictly greater than 12 cM and 220 SNP/cM.

```bash
workflows/07_ancibd/verify_and_extract_site_catalogs.sh \
  vcf_manifest.tsv work/site_catalogs
python3 workflows/07_ancibd/build_v62_resources.py \
  --v62 "$AADR_V62_SNP" --af-pattern "$AF_TABLE_PATTERN" \
  --vcf-sites-dir work/site_catalogs --output work/resources
workflows/07_ancibd/run_cross64.sh \
  vcf_manifest.tsv work/resources work/cross64
python3 workflows/07_ancibd/summarize_cross64.py \
  work/cross64/native_summary_0/ch_all.tsv \
  work/cross64/PAIR_DESIGN.tsv work/cross64
```

The VCF/BCF and ancIBD HDF5 intermediates are not distributed. Compact focal-pair
plotting summaries are included under `data/summary/ancibd/` and
`figures/_shared/ancibd_v62/`.
