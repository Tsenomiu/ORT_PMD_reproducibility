#!/usr/bin/env bash
# Path-neutral publication copy of the recovered ORT6–ORT27 sensitivity workflow.
set -euo pipefail

manifest=${1:?Usage: run_autosome_and_transversion_sensitivities.sh BAM_MANIFEST OUTPUT_DIR}
output_dir=${2:?Usage: run_autosome_and_transversion_sensitivities.sh BAM_MANIFEST OUTPUT_DIR}

: "${HG19_FASTA:?Set HG19_FASTA to the hg19 reference FASTA}"
: "${AADR_1240K_BED:?Set AADR_1240K_BED to the AADR v62.0 1240k BED}"
: "${AADR_V62_SNP:?Set AADR_V62_SNP to v62.0_1240k_public.snp}"
: "${READ2_PY:?Set READ2_PY to READv2 tag 2.1.0 READ2.py}"

SAMTOOLS=${SAMTOOLS:-samtools}
PILEUPCALLER=${PILEUPCALLER:-pileupCaller}
PLINK=${PLINK:-plink}
PYTHON=${PYTHON:-python3}
THREADS=${THREADS:-8}
MEMORY_MB=${MEMORY_MB:-20000}

for path in "$HG19_FASTA" "$AADR_1240K_BED" "$AADR_V62_SNP" "$READ2_PY"; do
  [[ -s "$path" ]] || { echo "Required file not found: $path" >&2; exit 1; }
done
[[ -s "$manifest" ]] || { echo "Manifest not found: $manifest" >&2; exit 1; }

read2_dir=$(cd "$(dirname "$READ2_PY")" && pwd -P)
READ2_PY="${read2_dir}/$(basename "$READ2_PY")"
manifest_dir=$(cd "$(dirname "$manifest")" && pwd -P)
manifest="${manifest_dir}/$(basename "$manifest")"
mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd -P)

samples=()
bams=()
while IFS=$'\t' read -r sample bam _rest; do
  [[ -z "$sample" || "$sample" == "sample_id" || "$sample" == \#* ]] && continue
  samples+=("$sample")
  bams+=("$bam")
done < "$manifest"

[[ ${#samples[@]} -eq 22 ]] || {
  echo "Expected 22 manifest rows ordered ORT6–ORT27." >&2
  exit 1
}
for i in "${!samples[@]}"; do
  expected="ORT$((i + 6))"
  [[ "${samples[$i]}" == "$expected" ]] || {
    echo "Manifest row $((i + 1)) must be $expected, not ${samples[$i]}." >&2
    exit 1
  }
  [[ -s "${bams[$i]}" ]] || { echo "BAM not found for $expected" >&2; exit 1; }
done
sample_names=$(IFS=,; printf '%s' "${samples[*]}")

pileup="$output_dir/cemetery.pileup.gz"
[[ ! -e "$pileup" ]] || { echo "Refusing to overwrite $pileup" >&2; exit 1; }

# Exact recovered pileup options: BAQ disabled, mapping/base quality 30,
# read-group tags ignored, and AADR-v62.0 1240k targets.
"$SAMTOOLS" mpileup -B -q30 -Q30 -R -l "$AADR_1240K_BED" \
  -f "$HG19_FASTA" "${bams[@]}" \
  | sed 's/^chr//' \
  | gzip -c > "$pileup"

# Exact recovered pseudo-haploid options. The original command set no seed.
gzip -dc "$pileup" \
  | "$PILEUPCALLER" --randomHaploid --sampleNames "$sample_names" \
      --samplePopName ORT -f "$AADR_V62_SNP" -p "$output_dir/cemetery_all"
gzip -dc "$pileup" \
  | "$PILEUPCALLER" --randomHaploid --skipTransitions \
      --sampleNames "$sample_names" --samplePopName ORT \
      -f "$AADR_V62_SNP" -p "$output_dir/cemetery_tv"

run_readv2() {
  local source_prefix=$1
  local label=$2
  local run_dir="$output_dir/read2_${label}"
  mkdir -p "$run_dir"
  "$PLINK" --bfile "$output_dir/$source_prefix" --chr 1-22 --make-bed \
    --allow-no-sex --threads "$THREADS" --memory "$MEMORY_MB" \
    --out "$run_dir/cemetery"
  (
    cd "$run_dir"
    "$PYTHON" "$READ2_PY" -i cemetery
  )
}

run_readv2 cemetery_all allsites
run_readv2 cemetery_tv transversions

echo "Autosome and transversion READv2 sensitivities complete: $output_dir"
