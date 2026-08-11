#!/usr/bin/env bash
# Publication wrapper reconstructed from the two recovered imputed-query scripts.
# The original 16-query and extra-8 merge order is retained.
set -euo pipefail
export LC_ALL=C

usage() {
  echo "Usage: $0 QUERY_MANIFEST AADR_1240K_BIM SCAFFOLD_PREFIX POSITION_TO_RSID OUTPUT_PREFIX WORK_DIR" >&2
  exit 2
}

[[ $# -eq 6 ]] || usage
manifest=$1
aadr_bim=$2
scaffold=$3
position_to_rsid=$4
output_prefix=$5
workdir=$6

PLINK=${PLINK:-plink}
BCFTOOLS=${BCFTOOLS:-bcftools}

for path in "$manifest" "$aadr_bim" "$position_to_rsid" \
            "${scaffold}.bed" "${scaffold}.bim" "${scaffold}.fam"; do
  [[ -f "$path" ]] || { echo "Missing input: $path" >&2; exit 2; }
done
[[ ! -e "$workdir" ]] || { echo "Work directory already exists: $workdir" >&2; exit 2; }
for suffix in bed bim fam; do
  [[ ! -e "${output_prefix}.${suffix}" ]] || { echo "Output already exists: ${output_prefix}.${suffix}" >&2; exit 2; }
done
mkdir -p "$workdir" "$(dirname "$output_prefix")"
awk 'BEGIN{OFS="\t"}{print $1,$4}' "$aadr_bim" > "$workdir/sites1240k.txt"

samples=()
vcfs=()
groups=()
while IFS=$'\t' read -r sample vcf _bam _rest; do
  [[ "$sample" == "sample_id" || -z "$sample" ]] && continue
  [[ -f "$vcf" ]] || { echo "Missing imputed VCF for $sample: $vcf" >&2; exit 2; }
  case "$sample" in
    ORT1[56]_fu_raw|ORT1[56]_fu_bamrefine5|ORT1[56]_fu_rescale5|ORT1[56]_fu_trim5|\
    ORT1[56]_nu_raw|ORT1[56]_nu_bamrefine5|ORT1[56]_nu_rescale5|ORT1[56]_nu_trim5)
      group=base16 ;;
    ORT1[56]_nu_bamrefine10|ORT1[56]_nu_rescale10|ORT1[56]_nu_rescaled|ORT1[56]_nu_trim10)
      group=extra8 ;;
    *) echo "Unexpected treatment state: $sample" >&2; exit 2 ;;
  esac
  samples+=("$sample")
  vcfs+=("$vcf")
  groups+=("$group")
done < "$manifest"

[[ ${#samples[@]} -eq 24 ]] || {
  echo "Expected 24 manifest rows; found ${#samples[@]}" >&2
  exit 2
}
expected=(
  ORT15_fu_raw ORT15_fu_bamrefine5 ORT15_fu_rescale5 ORT15_fu_trim5
  ORT15_nu_raw ORT15_nu_bamrefine5 ORT15_nu_bamrefine10 ORT15_nu_rescale5
  ORT15_nu_rescale10 ORT15_nu_rescaled ORT15_nu_trim5 ORT15_nu_trim10
  ORT16_fu_raw ORT16_fu_bamrefine5 ORT16_fu_rescale5 ORT16_fu_trim5
  ORT16_nu_raw ORT16_nu_bamrefine5 ORT16_nu_bamrefine10 ORT16_nu_rescale5
  ORT16_nu_rescale10 ORT16_nu_rescaled ORT16_nu_trim5 ORT16_nu_trim10
)
printf '%s\n' "${expected[@]}" | sort > "$workdir/expected_samples.txt"
printf '%s\n' "${samples[@]}" | sort > "$workdir/observed_samples.txt"
if ! diff -u "$workdir/expected_samples.txt" "$workdir/observed_samples.txt"; then
  echo "Manifest treatment states do not match the recovered 24-state design" >&2
  exit 2
fi

: > "$workdir/base16.inputs"
: > "$workdir/extra8.inputs"
for ((index=0; index<${#samples[@]}; index++)); do
  sample=${samples[$index]}
  vcf=${vcfs[$index]}
  group=${groups[$index]}
  prefix="$workdir/$sample"
  "$BCFTOOLS" view -T "$workdir/sites1240k.txt" "$vcf" -Oz \
    -o "${prefix}.sub.vcf.gz"
  "$PLINK" --vcf "${prefix}.sub.vcf.gz" --set-missing-var-ids @_# \
    --keep-allele-order --make-bed --allow-no-sex --out "$prefix"
  awk -v name="$sample" 'BEGIN{OFS="\t"}{$1=name;$2=name;print}' \
    "${prefix}.fam" > "${prefix}.renamed.fam"
  mv "${prefix}.renamed.fam" "${prefix}.fam"
  printf '%s\n' "$prefix" >> "$workdir/${group}.inputs"
done

merge_prefixes() {
  local input_list=$1
  local merged=$2
  local first=
  local merge_list="${merged}.merge-list"
  : > "$merge_list"
  while IFS= read -r prefix; do
    if [[ -z "$first" ]]; then
      first=$prefix
    else
      printf '%s.bed %s.bim %s.fam\n' "$prefix" "$prefix" "$prefix" >> "$merge_list"
    fi
  done < "$input_list"
  [[ -n "$first" ]] || { echo "Empty merge list: $input_list" >&2; exit 3; }
  set +e
  "$PLINK" --bfile "$first" --merge-list "$merge_list" --keep-allele-order \
    --make-bed --allow-no-sex --out "$merged"
  local status=$?
  set -e
  if [[ $status -ne 0 && -s "${merged}-merge.missnp" ]]; then
    local retry_list="${merged}.retry-list"
    local retry_first=
    : > "$retry_list"
    while IFS= read -r prefix; do
      local retry="${prefix}_aligned"
      "$PLINK" --bfile "$prefix" --exclude "${merged}-merge.missnp" \
        --keep-allele-order --make-bed --allow-no-sex --out "$retry"
      if [[ -z "$retry_first" ]]; then
        retry_first=$retry
      else
        printf '%s.bed %s.bim %s.fam\n' "$retry" "$retry" "$retry" >> "$retry_list"
      fi
    done < "$input_list"
    "$PLINK" --bfile "$retry_first" --merge-list "$retry_list" \
      --keep-allele-order --make-bed --allow-no-sex --out "$merged"
  elif [[ $status -ne 0 ]]; then
    echo "PLINK merge failed without a missnp list: $merged" >&2
    exit "$status"
  fi
}

harmonise_group() {
  local group=$1
  local expected_count=$2
  local merged="$workdir/${group}_merged"
  local scaffold_merge="$workdir/${group}_scaffold"
  local keep="$workdir/${group}.keep"
  local chrpos="$workdir/${group}_chrpos"
  local rsid="$workdir/${group}_rsid"
  merge_prefixes "$workdir/${group}.inputs" "$merged"

  set +e
  "$PLINK" --bfile "$merged" \
    --bmerge "${scaffold}.bed" "${scaffold}.bim" "${scaffold}.fam" \
    --keep-allele-order --make-bed --allow-no-sex --out "$scaffold_merge"
  local status=$?
  set -e
  if [[ $status -ne 0 && -s "${scaffold_merge}-merge.missnp" ]]; then
    "$PLINK" --bfile "$merged" --exclude "${scaffold_merge}-merge.missnp" \
      --keep-allele-order --make-bed --allow-no-sex --out "${merged}_aligned"
    "$PLINK" --bfile "$scaffold" --exclude "${scaffold_merge}-merge.missnp" \
      --keep-allele-order --make-bed --allow-no-sex --out "$workdir/scaffold_${group}_aligned"
    "$PLINK" --bfile "${merged}_aligned" \
      --bmerge "$workdir/scaffold_${group}_aligned.bed" \
               "$workdir/scaffold_${group}_aligned.bim" \
               "$workdir/scaffold_${group}_aligned.fam" \
      --keep-allele-order --make-bed --allow-no-sex --out "$scaffold_merge"
  elif [[ $status -ne 0 ]]; then
    echo "Scaffold merge failed without a missnp list: $group" >&2
    exit "$status"
  fi

  while IFS= read -r prefix; do
    sample=${prefix##*/}
    printf '%s %s\n' "$sample" "$sample"
  done < "$workdir/${group}.inputs" > "$keep"
  "$PLINK" --bfile "$scaffold_merge" --keep "$keep" \
    --extract <(awk '{print $2}' "${scaffold}.bim") --keep-allele-order \
    --make-bed --allow-no-sex --out "$chrpos"
  "$PLINK" --bfile "$chrpos" --update-name "$position_to_rsid" \
    --keep-allele-order --make-bed --allow-no-sex --out "$rsid"
  [[ $(wc -l < "${rsid}.fam") -eq "$expected_count" ]] || {
    echo "$group output sample count differs from $expected_count" >&2
    exit 3
  }
}

harmonise_group base16 16
harmonise_group extra8 8

set +e
"$PLINK" --bfile "$workdir/base16_rsid" \
  --bmerge "$workdir/extra8_rsid.bed" "$workdir/extra8_rsid.bim" \
           "$workdir/extra8_rsid.fam" \
  --keep-allele-order --make-bed --allow-no-sex --out "$output_prefix"
status=$?
set -e
if [[ $status -ne 0 && -s "${output_prefix}-merge.missnp" ]]; then
  "$PLINK" --bfile "$workdir/base16_rsid" --exclude "${output_prefix}-merge.missnp" \
    --keep-allele-order --make-bed --allow-no-sex --out "$workdir/base16_final"
  "$PLINK" --bfile "$workdir/extra8_rsid" --exclude "${output_prefix}-merge.missnp" \
    --keep-allele-order --make-bed --allow-no-sex --out "$workdir/extra8_final"
  "$PLINK" --bfile "$workdir/base16_final" \
    --bmerge "$workdir/extra8_final.bed" "$workdir/extra8_final.bim" \
             "$workdir/extra8_final.fam" \
    --keep-allele-order --make-bed --allow-no-sex --out "$output_prefix"
elif [[ $status -ne 0 ]]; then
  echo "Final 24-query merge failed without a missnp list" >&2
  exit "$status"
fi

[[ $(wc -l < "${output_prefix}.fam") -eq 24 ]] || {
  echo "Imputed output does not contain 24 samples" >&2
  exit 3
}
echo "Imputed query preparation complete: $output_prefix"
echo "Samples: 24; markers: $(wc -l < "${output_prefix}.bim")"
