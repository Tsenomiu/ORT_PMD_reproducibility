#!/usr/bin/env bash
# Run the ANGSD v0.933 male-X contamination analysis.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
bam=${1:?Usage: run_xchr_contamination.sh INPUT_BAM OUTPUT_PREFIX}
prefix=${2:?}
mkdir -p "$(dirname "$prefix")"

"$ANGSD" -i "$bam" -r X:5000000-154900000 \
  -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 20 \
  -ref "$HS37D5_FASTA" -out "$prefix"
"$ANGSD_CONTAMINATION" -a "$prefix.icnts.gz" -h "$HAPMAP_CHRX_SITES" -p "$THREADS" \
  > "$prefix.contamination.txt"
