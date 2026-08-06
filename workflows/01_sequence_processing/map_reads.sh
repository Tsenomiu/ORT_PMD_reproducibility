#!/usr/bin/env bash
# Map one single-end or paired-end read unit with BWA aln.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

mode=${1:?Usage: map_reads.sh SE|PE READ1 READ2_OR_DASH OUTPUT_BAM}
read1=${2:?}
read2=${3:?}
output=${4:?}
mkdir -p "$(dirname "$output")"
tmp="${output%.bam}.tmp"

case "$mode" in
  SE)
    "$BWA" aln -t "$THREADS" -l 16500 -n 0.01 -o 2 "$HS37D5_FASTA" "$read1" > "$tmp.sai"
    "$BWA" samse "$HS37D5_FASTA" "$tmp.sai" "$read1" \
      | "$SAMTOOLS" view -b -o "$output" -
    rm -f "$tmp.sai"
    ;;
  PE)
    [[ "$read2" != - ]] || { echo "PE mode needs READ2" >&2; exit 2; }
    "$BWA" aln -t "$THREADS" -l 16500 -n 0.01 -o 2 "$HS37D5_FASTA" "$read1" > "$tmp.1.sai"
    "$BWA" aln -t "$THREADS" -l 16500 -n 0.01 -o 2 "$HS37D5_FASTA" "$read2" > "$tmp.2.sai"
    "$BWA" sampe "$HS37D5_FASTA" "$tmp.1.sai" "$tmp.2.sai" "$read1" "$read2" \
      | "$SAMTOOLS" view -b -o "$output" -
    rm -f "$tmp.1.sai" "$tmp.2.sai"
    ;;
  *) echo "MODE must be SE or PE" >&2; exit 2 ;;
esac
