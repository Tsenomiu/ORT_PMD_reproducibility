#!/usr/bin/env bash
# Project the two ORT query representations onto shared AADR reference axes.
set -euo pipefail
source "${ORT_CONFIG:?Set ORT_CONFIG to a configured config.sh}"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
imputed=${1:?Usage: run_matched48_projection.sh IMP_PREFIX PHAP_PREFIX REF_PREFIX IID_POP KEEP_POPS OUTPUT_DIR}
phap=${2:?}; reference=${3:?}; iid_pop=${4:?}; keep_pops=${5:?}; outdir=${6:?}
mkdir -p "$outdir"; cd "$outdir"
[[ $(wc -l < "$imputed.fam") -eq 24 && $(wc -l < "$phap.fam") -eq 24 && $(wc -l < "$reference.fam") -eq 985 ]] || {
  echo "Expected 24 imputed, 24 pseudo-haploid, 985 reference individuals" >&2; exit 1;
}

"$PLINK2" --bfile "$reference" --chr 1-22 --set-all-var-ids @:# --make-bed --out aadr_ids
"$PLINK2" --bfile "$imputed" --chr 1-22 --set-all-var-ids @:# --make-bed --out imp_ids
"$PLINK2" --bfile "$phap" --chr 1-22 --set-all-var-ids @:# --make-bed --out phap_ids
"$PYTHON" "$script_dir/harmonize_markers.py" aadr_ids.bim imp_ids.bim phap_ids.bim .
[[ $(wc -l < common_sites_matched48.txt) -eq 1134702 ]] || { echo "Expected 1,134,702 aligned sites" >&2; exit 1; }

awk '{print $1,$2,"Projected_IMP","IMP__"$2}' imp_ids.fam > update_imp.txt
awk '{print $1,$2,"Projected_PHAP","PHAP__"$2}' phap_ids.fam > update_phap.txt
imp_flip=(); phap_flip=()
[[ -s flip_imp.txt ]] && imp_flip=(--flip flip_imp.txt)
[[ -s flip_phap.txt ]] && phap_flip=(--flip flip_phap.txt)
"$PLINK" --bfile aadr_ids --keep-allele-order --allow-no-sex --extract common_sites_matched48.txt --make-bed --out ref_c
"$PLINK" --bfile imp_ids --keep-allele-order --allow-no-sex "${imp_flip[@]}" --extract common_sites_matched48.txt --update-ids update_imp.txt --make-bed --out imp_c
"$PLINK" --bfile phap_ids --keep-allele-order --allow-no-sex "${phap_flip[@]}" --extract common_sites_matched48.txt --update-ids update_phap.txt --make-bed --out phap_c
printf 'imp_c\nphap_c\n' > merge.list
set +e; "$PLINK" --bfile ref_c --keep-allele-order --allow-no-sex --merge-list merge.list --make-bed --out merged; merge_status=$?; set -e
if [[ "$merge_status" -ne 0 && -s merged-merge.missnp ]]; then
  for prefix in ref_c imp_c phap_c; do "$PLINK" --bfile "$prefix" --exclude merged-merge.missnp --make-bed --out "${prefix}_x"; done
  printf 'imp_c_x\nphap_c_x\n' > merge2.list
  "$PLINK" --bfile ref_c_x --keep-allele-order --allow-no-sex --merge-list merge2.list --make-bed --out merged
elif [[ "$merge_status" -ne 0 ]]; then exit "$merge_status"; fi

# convertf treats missing PLINK phenotype as ignored; give every row a transient value.
awk '{$6=1; print}' OFS=' ' merged.fam > merged_pheno.fam
cat > convert.par <<EOF
genotypename:    merged.bed
snpname:         merged.bim
indivname:       merged_pheno.fam
outputformat:    EIGENSTRAT
genotypeoutname: matched48.geno
snpoutname:      matched48.snp
indivoutname:    matched48.ind
familynames:     NO
EOF
"$CONVERTF" -p convert.par > convertf.log 2>&1
"$PYTHON" "$script_dir/label_eigenstrat.py" matched48.ind "$iid_pop" matched48.labelled.ind
cat > smartpca.par <<EOF
genotypename:    matched48.geno
snpname:         matched48.snp
indivname:       matched48.labelled.ind
evecoutname:     pca_aadr_matched48.evec
evaloutname:     pca_aadr_matched48.eval
poplistname:     $keep_pops
numoutevec:      10
numthreads:      $THREADS
lsqproject:      YES
numoutlieriter:  0
altnormstyle:    NO
EOF
"$SMARTPCA" -p smartpca.par > smartpca.log 2>&1
grep -q 'total number of snps killed in pass: 20641  used: 1113348' smartpca.log || {
  echo "smartpca internal-site count differs from the expected study value" >&2; exit 1;
}
