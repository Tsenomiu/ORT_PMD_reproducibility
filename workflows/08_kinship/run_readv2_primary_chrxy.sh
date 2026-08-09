#!/usr/bin/env bash
# Reconstruct the calling logic for the primary ORT1--ORT27 READv2 analysis.
# The 1240k BED/SNP resources must include chromosomes 1--22, X, and Y.
# Inputs are merged full-UDG a--c BAMs after 6-bp terminal soft-clipping.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
manifest=${1:?Usage: run_readv2_primary_chrxy.sh BAM_MANIFEST OUTPUT_DIR}
outdir=${2:?}
manifest="$(cd "$(dirname "$manifest")" && pwd -P)/$(basename "$manifest")"
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd -P)
samples=(); bams=()
while IFS=$'\t' read -r sample bam _rest; do
  [[ -z "$sample" || "$sample" == sample || "$sample" == \#* ]] && continue
  samples[${#samples[@]}]="$sample"
  bams[${#bams[@]}]="$bam"
done < "$manifest"
[[ ${#samples[@]} -eq 27 && ${#bams[@]} -eq 27 ]] || { echo "Expected ORT1--ORT27 (27 rows)" >&2; exit 1; }
for i in "${!samples[@]}"; do
  [[ "${samples[$i]}" == "ORT$((i+1))" && -s "${bams[$i]}" ]] || { echo "Manifest must be ordered ORT1--ORT27" >&2; exit 1; }
done
names=$(IFS=,; printf '%s' "${samples[*]}")

"$SAMTOOLS" mpileup -B -q 30 -Q 30 -R -l "$AADR_1240K_BED" -f "$HG19_FASTA" "${bams[@]}" \
  | sed 's/^chr//' | gzip > "$outdir/cemetery_chrxy.pileup.gz"
gzip -dc "$outdir/cemetery_chrxy.pileup.gz" \
  | "$PILEUPCALLER" --randomHaploid --sampleNames "$names" --samplePopName ORT \
      -f "$AADR_V62_SNP" -p "$outdir/cemetery_1240k_chrxy"
(cd "$outdir" && "$PYTHON" "$READ2_PY" -i "cemetery_1240k_chrxy")
