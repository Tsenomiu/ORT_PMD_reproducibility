# 01 Sequence processing

This stage removes adapters, maps reads to hs37d5, merges sequencing runs, and
removes PCR duplicates.

- `run_adapterremoval.sh` uses AdapterRemoval v2.2.2 with trimming of ambiguous
  and low-quality ends, a 30-bp minimum length, minimum quality 25, and read-pair
  collapsing.
- `map_reads.sh` uses BWA aln (`-l 16500 -n 0.01 -o 2`) and selects `samse` or
  `sampe` according to read type.
- `merge_markdup.sh` merges BAMs and applies the appropriate samtools duplicate
  workflow.

Three read sets can be assembled: collapsed reads only; collapsed plus
`collapsed.truncated`; or both collapsed classes plus unmerged paired reads. The
reported imputation and IBD analyses use collapsed reads only.

```bash
export ORT_CONFIG="$PWD/workflows/config.sh"
workflows/01_sequence_processing/run_adapterremoval.sh fastq_manifest.tsv work/trimmed
workflows/01_sequence_processing/map_reads.sh SE work/trimmed/x.collapsed.gz - work/mapped/x.bam
workflows/01_sequence_processing/merge_markdup.sh mapped_bam_manifest.tsv work/deduplicated
```

Required FASTQ files are obtained from ENA and are not stored here.
