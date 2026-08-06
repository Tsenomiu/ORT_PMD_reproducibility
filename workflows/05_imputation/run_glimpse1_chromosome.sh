#!/usr/bin/env bash
# Run the GLIMPSE v1 workflow for one chromosome and treatment.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

sample=${1:?Usage: run_glimpse1_chromosome.sh SAMPLE BAM CHROM OUTPUT_ROOT}
bam=${2:?}
chrom=${3:?}
root=${4:?}
subst_chr() { printf '%s\n' "${1//\{CHROM\}/$chrom}"; }
ref_vcf=$(subst_chr "$PHASE3_REF_PATTERN")
sites_vcf=$(subst_chr "$PHASE3_SITES_VCF_PATTERN")
sites_tsv=$(subst_chr "$PHASE3_SITES_TSV_PATTERN")
gmap=$(subst_chr "$GENETIC_MAP_PATTERN")
out="$root/$sample/chr$chrom"
mkdir -p "$out/imputed"

for path in "$HS37D5_FASTA" "$bam" "$ref_vcf" "$sites_vcf" "$sites_tsv" "$gmap"; do
  [[ -s "$path" ]] || { echo "Missing input: $path" >&2; exit 1; }
done

glvcf="$out/$sample.chr$chrom.gl.vcf.gz"
chunks="$out/chunks.chr$chrom.txt"
"$BCFTOOLS" mpileup -f "$HS37D5_FASTA" --ignore-RG -I -E -a FORMAT/DP \
  -T "$ref_vcf" -r "$chrom" -q 30 -Q 20 "$bam" -Ou --threads "$THREADS" \
  | "$BCFTOOLS" call -Aim -C alleles -T "$sites_tsv" -Oz -o "$glvcf" --threads "$THREADS"
"$TABIX" -f "$glvcf"

"$GLIMPSE1_CHUNK" --input "$glvcf" --region "$chrom" \
  --window-size 2000000 --buffer-size 200000 --output "$chunks"

list="$out/list.chr$chrom.txt"
: > "$list"
while read -r chunk_id _ input_region output_region _; do
  id=$(printf '%02d' "$chunk_id")
  bcf="$out/imputed/$sample.chr$chrom.$id.bcf"
  "$GLIMPSE1_PHASE" --input "$glvcf" --reference "$ref_vcf" --map "$gmap" \
    --input-region "$input_region" --output-region "$output_region" \
    --output "$bcf" --thread 1
  "$BCFTOOLS" index -f "$bcf"
  printf '%s\n' "$bcf" >> "$list"
done < "$chunks"

ligated="$out/$sample.chr$chrom.ligated.bcf"
phased="$out/$sample.chr$chrom.phased.bcf"
"$GLIMPSE1_LIGATE" --input "$list" --output "$ligated"
"$BCFTOOLS" index -f "$ligated"
"$GLIMPSE1_SAMPLE" --input "$ligated" --solve --output "$phased" --thread 1
"$BCFTOOLS" index -f "$phased"
