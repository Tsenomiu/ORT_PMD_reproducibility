#!/usr/bin/env bash
# Run ancIBD for the 64 ordered treatment comparisons.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
manifest=${1:?Usage: run_cross64.sh VCF_MANIFEST RESOURCE_DIR OUTPUT_DIR}
resources=${2:?}
outdir=${3:?}
mkdir -p "$outdir/chromosomes" "$outdir/logs"

version=$("$PYTHON" -m pip show ancIBD 2>/dev/null | awk '$1=="Version:"{print $2}')
[[ "$version" == 0.7 ]] || { echo "Expected ancIBD 0.7, found ${version:-none}" >&2; exit 1; }
cp "$manifest" "$outdir/SAMPLE_MANIFEST.tsv"
awk -F '\t' 'NR>1 && $1!="" && $1!~/^#/{print $1}' "$manifest" > "$outdir/iids.txt"
[[ $(wc -l < "$outdir/iids.txt") -eq 16 ]] || { echo "Expected 16 samples" >&2; exit 1; }
mapfile -t ids15 < <(awk -F '\t' 'NR>1 && $1~/^ORT15_/{print $1}' "$manifest")
mapfile -t ids16 < <(awk -F '\t' 'NR>1 && $1~/^ORT16_/{print $1}' "$manifest")
: > "$outdir/pairs_64.txt"; printf 'iid1\tiid2\n' > "$outdir/PAIR_DESIGN.tsv"
for iid1 in "${ids15[@]}"; do for iid2 in "${ids16[@]}"; do
  printf '%s %s\n' "$iid1" "$iid2" >> "$outdir/pairs_64.txt"
  printf '%s\t%s\n' "$iid1" "$iid2" >> "$outdir/PAIR_DESIGN.tsv"
done; done
[[ $(wc -l < "$outdir/pairs_64.txt") -eq 64 ]] || exit 1

map="$resources/map/v62.autosomes.Morgan.snp"
for chrom in $(seq 1 22); do
  chrdir="$outdir/chromosomes/ch$chrom"; tmp="$chrdir/tmp"
  mkdir -p "$tmp"; : > "$tmp/bcfs.list"
  marker="$resources/markers/v62_targets_ch${chrom}.tsv"
  af="$resources/afs/v62_1000G_AF_ch${chrom}.tsv"
  while IFS=$'\t' read -r iid vcf; do
    [[ -z "$iid" || "$iid" == iid || "$iid" == \#* ]] && continue
    raw="$tmp/$iid.raw.bcf"; renamed="$tmp/$iid.bcf"
    "$BCFTOOLS" view -R "$marker" -Ob -o "$raw" "$vcf"
    printf '%s\n' "$iid" > "$tmp/$iid.name"
    "$BCFTOOLS" reheader -s "$tmp/$iid.name" -o "$renamed" "$raw"
    "$BCFTOOLS" index -f "$renamed"; printf '%s\n' "$renamed" >> "$tmp/bcfs.list"
  done < "$manifest"
  cohort="$chrdir/cohort.ch$chrom.bcf"
  "$BCFTOOLS" merge -m none -Ob -l "$tmp/bcfs.list" -o "$cohort"
  "$BCFTOOLS" index -f "$cohort"
  [[ $("$BCFTOOLS" index -n "$cohort") -eq $(wc -l < "$marker") ]] || { echo "chr$chrom site mismatch" >&2; exit 1; }
  "$ANCIBD_RUN" --vcf "$cohort" --ch "$chrom" --marker_path "$marker" \
    --map_path "$map" --af_path "$af" --af_column variants/AF_ALL \
    --iid "$outdir/iids.txt" --pair "$outdir/pairs_64.txt" \
    --out "$outdir" --prefix p --min 8 --ibd-in 1 --ibd-out 10 -v \
    > "$outdir/logs/ch$chrom.log" 2>&1
done

mkdir -p "$outdir/native_summary_0" "$outdir/native_summary_220"
"$ANCIBD_SUMMARY" --tsv "$outdir/p.ch" --ch 1-22 --bin 8,12,16,20 --snp_cm 0 \
  --out "$outdir/native_summary_0" > "$outdir/logs/summary_0.log" 2>&1
"$ANCIBD_SUMMARY" --tsv "$outdir/p.ch" --ch 1-22 --bin 8,12,16,20 --snp_cm 220 \
  --out "$outdir/native_summary_220" > "$outdir/logs/summary_220.log" 2>&1
