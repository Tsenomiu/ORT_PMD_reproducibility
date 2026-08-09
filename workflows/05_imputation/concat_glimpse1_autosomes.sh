#!/usr/bin/env bash
# Concatenate the 22 normalised per-chromosome GLIMPSE outputs for one dataset.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

sample=${1:?Usage: concat_glimpse1_autosomes.sh DATASET OUTPUT_ROOT}
root=${2:?}
dataset_root="$root/$sample"
list="$dataset_root/$sample.autosomes.list.txt"
output="$dataset_root/$sample.imputed.allchr.vcf.gz"

: > "$list"
for chrom in $(seq 1 22); do
  vcf="$dataset_root/chr$chrom/$sample.chr$chrom.phased.withGP.normalised.vcf.gz"
  [[ -s "$vcf" && -s "$vcf.tbi" ]] || {
    echo "Missing phased chromosome input or index: $vcf" >&2
    exit 1
  }
  printf '%s\n' "$vcf" >> "$list"
done

"$BCFTOOLS" concat -f "$list" -Oz -o "$output" --threads "$THREADS"
"$TABIX" -f -p vcf "$output"

chromosomes=$("$BCFTOOLS" index -s "$output" | cut -f1 | tr '\n' ' ')
expected=$(seq 1 22 | tr '\n' ' ')
[[ "$chromosomes" == "$expected" ]] || {
  echo "Autosome index check failed: $chromosomes" >&2
  exit 1
}

printf 'wrote %s\n' "$output"
