#!/usr/bin/env bash
# Run standalone schmutzi 1.5.7 mtCont with asymmetric damage profiles.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
bam=${1:?Usage: run_mtdna_contamination.sh INPUT_NONUDG_BAM OUTPUT_DIR SAMPLE}
outdir=${2:?}
sample=${3:?}
mkdir -p "$outdir"

"$BAM2PROF" -length 5 -double \
  -5p "$outdir/$sample.5p.prof" -3p "$outdir/$sample.3p.prof" "$bam"
"$ENDOCALLER" -seq "$outdir/$sample.endo.fa" -log "$outdir/$sample.endo.log" \
  -name MT -qual 0 -logindel 50 -single -l 16569 "$RCRS_FASTA" "$bam"
"$MTCONT" -deam5p "$outdir/$sample.5p.prof" -deam3p "$outdir/$sample.3p.prof" \
  -o "$outdir/$sample.mtcont" "$outdir/$sample.endo.log" "$RCRS_FASTA" "$bam" \
  "$SCHMUTZI_FREQ_DIR"/*.freq
"$CONTOUT2CONTEST" "$outdir/$sample.mtcont" > "$outdir/$sample.contamination.txt"
