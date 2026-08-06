#!/usr/bin/env bash
# Call haploid mitochondrial variants with bcftools 1.16 at MAPQ/BQ 30/30.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
manifest=${1:?Usage: run_haplogroup.sh BAM_MANIFEST SAMPLE OUTPUT_DIR}
sample=${2:?}
outdir=${3:?}
mkdir -p "$outdir/fragments"
list="$outdir/$sample.fragments.list"; : > "$list"
index=0
while IFS=$'\t' read -r manifest_sample bam; do
  [[ -z "$manifest_sample" || "$manifest_sample" == "sample" || "$manifest_sample" == \#* ]] && continue
  [[ "$manifest_sample" == "$sample" && -s "$bam" ]] || continue
  index=$((index + 1)); fragment="$outdir/fragments/$sample.$index.MT_Y.bam"
  "$SAMTOOLS" view -b -o "$fragment" "$bam" MT Y
  printf '%s\n' "$fragment" >> "$list"
done < "$manifest"
[[ "$index" -eq 3 ]] || { echo "Expected three full-UDG sublibrary BAMs for $sample" >&2; exit 1; }
merged="$outdir/$sample.MT_Y.merged.bam"
"$SAMTOOLS" merge -f -b "$list" "$merged"
"$SAMTOOLS" index "$merged"
unlabelled="$outdir/$sample.mt.unlabelled.vcf.gz"
"$BCFTOOLS" mpileup -r MT -f "$HS37D5_FASTA" -B -q 30 -Q 30 "$merged" -Ou \
  | "$BCFTOOLS" call --ploidy 1 -m -Oz -o "$unlabelled"
printf '%s\n' "$sample" > "$outdir/$sample.sample_name.txt"
"$BCFTOOLS" reheader -s "$outdir/$sample.sample_name.txt" -o "$outdir/$sample.mt.vcf.gz" "$unlabelled"
"$TABIX" -f -p vcf "$outdir/$sample.mt.vcf.gz"
"$HAPLOGREP3" classify --in "$outdir/$sample.mt.vcf.gz" \
  --out "$outdir/$sample.haplogrep.txt" --tree phylotree-fu-rcrs@1.3
