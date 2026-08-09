#!/usr/bin/env bash
# Run KING-robust and IBS0 calculations from the ORT15 and ORT16 VCFs.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
vcf15=${1:?Usage: run_king_ibs0.sh ORT15_VCF ORT16_VCF OUTPUT_DIR}
vcf16=${2:?}
outdir=${3:?}
mkdir -p "$outdir"

"$BCFTOOLS" merge --force-samples -m none "$vcf15" "$vcf16" -Ou \
  | "$BCFTOOLS" view -i 'INFO/RAF[0]>=0.1 && INFO/RAF[0]<=0.9' -Ou \
  | "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/RAF[\t%GT\t%GP]\n' \
  > "$outdir/common_raf01_09.tsv"
"$PYTHON" "$script_dir/summarize_king_ibs0.py" "$outdir/common_raf01_09.tsv" \
  "$outdir/king_ibs0_genomewide.tsv"

"$BCFTOOLS" merge --force-samples -m none -r 1 "$vcf15" "$vcf16" -Ou \
  | "$BCFTOOLS" view -i 'INFO/RAF[0]>=0.1 && INFO/RAF[0]<=0.9' -Ou \
  | "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/RAF[\t%GT\t%GP]\n' \
  > "$outdir/common_raf01_09_chr1.tsv"
"$PYTHON" "$script_dir/summarize_king_ibs0.py" --min-max-gp 0.99 \
  "$outdir/common_raf01_09_chr1.tsv" "$outdir/king_ibs0_chr1_gp99.tsv"
