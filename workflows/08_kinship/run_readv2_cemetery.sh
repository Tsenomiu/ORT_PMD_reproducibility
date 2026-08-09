#!/usr/bin/env bash
# Reconstruct the later ORT6--ORT27 autosome-only READv2 sensitivity run.
# Inputs are merged full-UDG a--c BAMs after 6-bp terminal soft-clipping.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
manifest=${1:?Usage: run_readv2_cemetery.sh BAM_MANIFEST OUTPUT_DIR}
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
[[ ${#samples[@]} -eq 22 && ${#bams[@]} -eq 22 ]] || { echo "Expected ORT6--ORT27 (22 rows)" >&2; exit 1; }
for i in "${!samples[@]}"; do
  [[ "${samples[$i]}" == "ORT$((i+6))" && -s "${bams[$i]}" ]] || { echo "Manifest must be ordered ORT6--ORT27" >&2; exit 1; }
done
names=$(IFS=,; printf '%s' "${samples[*]}")

"$SAMTOOLS" mpileup -B -q 30 -Q 30 -R -l "$AADR_1240K_BED" -f "$HG19_FASTA" "${bams[@]}" \
  | sed 's/^chr//' | gzip > "$outdir/cemetery.pileup.gz"
gzip -dc "$outdir/cemetery.pileup.gz" \
  | "$PILEUPCALLER" --randomHaploid --sampleNames "$names" --samplePopName ORT \
      -f "$AADR_V62_SNP" -p "$outdir/cemetery_1240k"
"$PLINK" --bfile "$outdir/cemetery_1240k" --chr 1-22 --make-bed --allow-no-sex \
  --threads "$THREADS" --out "$outdir/cemetery_1240k_autosomes"
(cd "$outdir" && "$PYTHON" "$READ2_PY" -i "cemetery_1240k_autosomes")
