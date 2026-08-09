#!/usr/bin/env bash
# Count reference and alternate reads with the study filters.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
bam=${1:?Usage: run_alt_counts.sh BAM SITES_VCF SITES_ALLELES_TSV OUTPUT_TSV}
sites_vcf=${2:?}
alleles_tsv=${3:?}
output=${4:?}
mkdir -p "$(dirname "$output")"
{
  printf 'CHROM\tPOS\tREF\tALT\tDP\tAD_REF\tAD_ALT\n'
  "$BCFTOOLS" mpileup -f "$HS37D5_FASTA" --ignore-RG -I -E \
    -a FORMAT/DP,FORMAT/AD -T "$sites_vcf" -q 30 -Q 20 -d 10000 "$bam" -Ou \
    | "$BCFTOOLS" call -Aim -C alleles -T "$alleles_tsv" -Ou \
    | "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD{0}\t%AD{1}]\n'
} > "$output"
