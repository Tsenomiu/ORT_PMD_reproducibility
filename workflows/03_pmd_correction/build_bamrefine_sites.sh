#!/usr/bin/env bash
# Create the six-column SNP list required by bamRefine.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
vcf=${1:?Usage: build_bamrefine_sites.sh PHASE3_VCF OUTPUT_SNP_LIST}
output=${2:?}
mkdir -p "$(dirname "$output")"
"$BCFTOOLS" view -m2 -M2 -v snps "$vcf" -Ou \
  | "$BCFTOOLS" query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' \
  | awk 'BEGIN{OFS="\t"}{id=($1=="."?$2":"$3":"$4":"$5:$1); print id,$2,0,$3,$4,$5}' \
  > "$output"
