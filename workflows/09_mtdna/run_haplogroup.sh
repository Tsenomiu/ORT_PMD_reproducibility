#!/usr/bin/env bash
# Call haploid mitochondrial variants from one canonical collapsed-only BAM.
set -euo pipefail

source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"

manifest=${1:?Usage: run_haplogroup.sh BAM_MANIFEST SAMPLE OUTPUT_DIR}
sample=${2:?}
outdir=${3:?}

mkdir -p "$outdir/bams" "$outdir/vcf" "$outdir/haplogrep" "$outdir/provenance"

# The manifest validator permits the verified combined a--c collapsed-only state
# and verifies its known SHA-256, rejecting renamed or nested read-set supersets.
# In particular, .1 is the verified canonical input; .2 and .3 are not inputs.
bam=$(
  "$PYTHON" "$(dirname "${BASH_SOURCE[0]}")/validate_manifest.py" \
    "$manifest" --sample "$sample" --require-files
)
"$SAMTOOLS" quickcheck -v "$bam"

mt_bam="$outdir/bams/$sample.fullUDG_collapsed_only.MT.bam"
"$SAMTOOLS" view -@ "$THREADS" -bh -o "$mt_bam" "$bam" MT
"$SAMTOOLS" index -@ "$THREADS" "$mt_bam"
"$SAMTOOLS" quickcheck -v "$mt_bam"

"$SAMTOOLS" --version > "$outdir/provenance/samtools.version.txt"
"$BCFTOOLS" --version > "$outdir/provenance/bcftools.version.txt"
printf 'sample\tread_set_state\tinput_scope\tinput_bam\n%s\t%s\t%s\t%s\n' \
  "$sample" full_udg_collapsed_only combined_libraries_a_b_c "$bam" \
  > "$outdir/provenance/$sample.input.tsv"

# With bcftools 1.16, omitting -d retains the documented default maximum input
# depth of 250 reads per file. No minimum depth, allele-fraction or strand-balance
# filter is added here.
raw_vcf="$outdir/vcf/$sample.mtDNA.raw_sample_name.vcf.gz"
final_vcf="$outdir/vcf/$sample.mtDNA.collapsed_only.vcf.gz"
"$BCFTOOLS" mpileup -r MT -f "$HS37D5_FASTA" -B -q 30 -Q 30 "$mt_bam" -Ou \
  | "$BCFTOOLS" call --ploidy 1 -m -Oz -o "$raw_vcf"
"$BCFTOOLS" index -f -t "$raw_vcf"
printf '%s\n' "$sample" > "$outdir/provenance/$sample.sample_name.txt"
"$BCFTOOLS" reheader -s "$outdir/provenance/$sample.sample_name.txt" \
  -o "$final_vcf" "$raw_vcf"
"$BCFTOOLS" index -f -t "$final_vcf"

"$HAPLOGREP3" classify --in "$final_vcf" \
  --out "$outdir/haplogrep/$sample.haplogrep.txt" \
  --tree phylotree-fu-rcrs@1.3
