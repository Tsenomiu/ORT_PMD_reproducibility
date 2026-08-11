#!/usr/bin/env bash
# Path-neutral publication copy of the recovered 24-query pseudo-haploid command.
set -euo pipefail
export LC_ALL=C

usage() {
  echo "Usage: $0 QUERY_MANIFEST HS37D5_FASTA AADR_1240K_BED AADR_1240K_SNP OUTPUT_PREFIX WORK_DIR" >&2
  exit 2
}

[[ $# -eq 6 ]] || usage
manifest=$1
reference_fasta=$2
sites_bed=$3
sites_snp=$4
output_prefix=$5
workdir=$6

SAMTOOLS=${SAMTOOLS:-samtools}
PILEUPCALLER=${PILEUPCALLER:-pileupCaller}

for path in "$manifest" "$reference_fasta" "$sites_bed" "$sites_snp"; do
  [[ -f "$path" ]] || { echo "Missing input: $path" >&2; exit 2; }
done
[[ ! -e "$workdir" ]] || { echo "Work directory already exists: $workdir" >&2; exit 2; }
for suffix in bed bim fam; do
  [[ ! -e "${output_prefix}.${suffix}" ]] || { echo "Output already exists: ${output_prefix}.${suffix}" >&2; exit 2; }
done
mkdir -p "$workdir" "$(dirname "$output_prefix")"

samples=()
bams=()
while IFS=$'\t' read -r sample _imputed_vcf bam _rest; do
  [[ "$sample" == "sample_id" || -z "$sample" ]] && continue
  [[ -f "$bam" ]] || { echo "Missing BAM for $sample: $bam" >&2; exit 2; }
  samples+=("$sample")
  bams+=("$bam")
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

for bam in "${bams[@]}"; do
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.*}.bai" ]]; then
    "$SAMTOOLS" index "$bam"
  fi
done

sample_names=$(IFS=,; echo "${samples[*]}")
"$SAMTOOLS" mpileup -B -q30 -Q30 -R -l "$sites_bed" -f "$reference_fasta" \
  "${bams[@]}" 2> "$workdir/mpileup24.err" \
  | "$PILEUPCALLER" --randomHaploid --sampleNames "$sample_names" \
      --samplePopName ORT -f "$sites_snp" -p "$output_prefix" \
      2> "$workdir/pileupcaller24.err"

[[ $(wc -l < "${output_prefix}.fam") -eq 24 ]] || {
  echo "Pseudo-haploid output does not contain 24 samples" >&2
  exit 3
}
echo "Pseudo-haploid query preparation complete: $output_prefix"
echo "Samples: 24; markers: $(wc -l < "${output_prefix}.bim")"
