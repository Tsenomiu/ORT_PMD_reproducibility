#!/usr/bin/env bash
# Verify matching CHROM/POS/REF/ALT catalogs across 16 v62-restricted VCFs.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
manifest=${1:?Usage: verify_and_extract_site_catalogs.sh VCF_MANIFEST OUTPUT_DIR}
outdir=${2:?}
mkdir -p "$outdir/catalogs"
regions="$outdir/v62_autosomal_regions.tsv"
awk 'BEGIN{OFS="\t"} $2>=1 && $2<=22{print $2,$4}' "$AADR_V62_SNP" > "$regions"

[[ $(awk -F '\t' 'NR>1 && $1!="" && $1!~/^#/{n++}END{print n+0}' "$manifest") -eq 16 ]] || { echo "Expected 16 VCFs" >&2; exit 1; }
first=1
while IFS=$'\t' read -r iid vcf; do
  [[ -z "$iid" || "$iid" == iid || "$iid" == \#* ]] && continue
  [[ -s "$vcf" ]] || { echo "Missing $vcf" >&2; exit 1; }
  sample_dir="$outdir/catalogs/$iid"
  mkdir -p "$sample_dir"
  for chrom in $(seq 1 22); do
    "$BCFTOOLS" query -r "$chrom" -R "$regions" -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" \
      > "$sample_dir/vcf_sites_ch${chrom}.tsv"
    if [[ "$first" -eq 1 ]]; then
      cp "$sample_dir/vcf_sites_ch${chrom}.tsv" "$outdir/vcf_sites_ch${chrom}.tsv"
    else
      cmp -s "$outdir/vcf_sites_ch${chrom}.tsv" "$sample_dir/vcf_sites_ch${chrom}.tsv" || {
        echo "Site catalog mismatch for $iid chromosome $chrom" >&2; exit 1;
      }
    fi
  done
  first=0
done < "$manifest"

rows=$(awk 'END{print NR}' "$outdir"/vcf_sites_ch*.tsv | awk '{s+=$1}END{print s}')
[[ "$rows" -gt 1000000 ]] || { echo "Implausibly few common sites: $rows" >&2; exit 1; }
printf 'status=PASS\ninputs=16\ncomparison_key=CHROM_POS_REF_ALT\nrows=%s\n' "$rows" > "$outdir/PASS.txt"
# Retain the per-input catalogs with the run outputs for site-intersection checks.
