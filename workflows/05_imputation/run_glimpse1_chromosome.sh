#!/usr/bin/env bash
# Run the GLIMPSE v1 workflow for one chromosome and treatment.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

sample=${1:?Usage: run_glimpse1_chromosome.sh DATASET INDIVIDUAL_ID BAM CHROM OUTPUT_ROOT}
individual=${2:?}
bam=${3:?}
chrom=${4:?}
root=${5:?}
subst_chr() { printf '%s\n' "${1//\{CHROM\}/$chrom}"; }
ref_vcf=$(subst_chr "$PHASE3_REF_PATTERN")
sites_tsv=$(subst_chr "$PHASE3_SITES_TSV_PATTERN")
gmap=$(subst_chr "$GENETIC_MAP_PATTERN")
out="$root/$sample/chr$chrom"
mkdir -p "$out/imputed"

for path in "$HS37D5_FASTA" "$bam" "$ref_vcf" "$sites_tsv" "$gmap"; do
  [[ -s "$path" ]] || { echo "Missing input: $path" >&2; exit 1; }
done

glvcf="$out/$sample.chr$chrom.gl.vcf.gz"
chunks="$out/chunks.chr$chrom.txt"
"$BCFTOOLS" mpileup -f "$HS37D5_FASTA" --ignore-RG -I -E -a FORMAT/DP \
  -T "$ref_vcf" -r "$chrom" -q 30 -Q 20 -d 250 "$bam" -Ou --threads "$THREADS" \
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
    --seed 15052011 --output "$bcf" --thread 1
  "$BCFTOOLS" index -f "$bcf"
  printf '%s\n' "$bcf" >> "$list"
done < "$chunks"

ligated="$out/$sample.chr$chrom.ligated.bcf"
phased="$out/$sample.chr$chrom.phased.bcf"
"$GLIMPSE1_LIGATE" --input "$list" --output "$ligated"
"$BCFTOOLS" index -f "$ligated"
"$GLIMPSE1_SAMPLE" --input "$ligated" --solve --output "$phased" --thread 1
"$BCFTOOLS" index -f "$phased"

# GLIMPSE_sample --solve writes phased GT but not GP. Restore GP from the
# position- and sample-matched ligated output, then give every treatment from an
# individual the same sample identifier required by the concordance workflow.
withgp="$out/$sample.chr$chrom.phased.withGP.vcf.gz"
normalised="$out/$sample.chr$chrom.phased.withGP.normalised.vcf.gz"
sample_names="$out/$sample.chr$chrom.sample_names.txt"
"$BCFTOOLS" annotate -a "$ligated" -c +FORMAT/GP -Oz -o "$withgp" \
  --threads "$THREADS" "$phased"
"$TABIX" -f -p vcf "$withgp"
printf '%s\n' "$individual" > "$sample_names"
"$BCFTOOLS" reheader -s "$sample_names" -o "$normalised" "$withgp"
"$TABIX" -f -p vcf "$normalised"

observed=$("$BCFTOOLS" query -l "$normalised")
[[ "$observed" == "$individual" ]] || {
  echo "Unexpected sample ID in $normalised: $observed" >&2
  exit 1
}

# Validate every record, rather than relying on the generic GT description in
# the GLIMPSE header. --solve must yield a diploid phased GT, and the copied GP
# field must contain three posterior probabilities for concordance.
if ! "$BCFTOOLS" query -f '[%GT\t%GP\n]' "$normalised" | awk -F '\t' '
  BEGIN { records = 0; invalid = 0 }
  {
    records++
    if ($1 !~ /^[0-9]+[|][0-9]+$/) invalid++
    n_gp = split($2, gp, ",")
    if (n_gp != 3 || gp[1] == "." || gp[2] == "." || gp[3] == ".") invalid++
  }
  END { exit(records == 0 || invalid != 0) }
'; then
  echo "Phased-GT/GP validation failed: $normalised" >&2
  exit 1
fi
