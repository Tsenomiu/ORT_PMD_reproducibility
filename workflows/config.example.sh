#!/usr/bin/env bash
# Copy to config.sh and edit locally. config.sh is ignored by Git.
# Use absolute paths on the machine that will run the analysis.

export THREADS="${THREADS:-8}"
export PYTHON="${PYTHON:-python3}"
export WORK_ROOT="${WORK_ROOT:-/path/to/analysis_work}"
export PRIVATE_DATA_ROOT="${PRIVATE_DATA_ROOT:-/path/to/restricted_input_data}"

# References and site panels
export HS37D5_FASTA="${HS37D5_FASTA:-/path/to/hs37d5.fa}"
export HG19_FASTA="${HG19_FASTA:-/path/to/hg19.fa}"
export RCRS_FASTA="${RCRS_FASTA:-/path/to/rCRS_NC_012920.fa}"
export HAPMAP_CHRX_SITES="${HAPMAP_CHRX_SITES:-/path/to/HapMapChrX.gz}"
export SCHMUTZI_FREQ_DIR="${SCHMUTZI_FREQ_DIR:-/path/to/schmutzi/alleleFreqMT/eurasian/freqs}"
export BAMREFINE_SNP_LIST="${BAMREFINE_SNP_LIST:-/path/to/1000G_biallelic_GRCh37.bamrefine.snp}"
export AADR_V62_SNP="${AADR_V62_SNP:-/path/to/v62.0_1240k_public.snp}"
export AADR_1240K_BED="${AADR_1240K_BED:-/path/to/v62_1240k_sites.bed}"
export AADR_REFERENCE_PREFIX="${AADR_REFERENCE_PREFIX:-/path/to/aadr_v62_presentday_hgdp_sgdp}"
export PCA_REFERENCE_POPULATIONS="${PCA_REFERENCE_POPULATIONS:-/path/to/presentday_hgdp_sgdp.populations.txt}"

# `{CHROM}` is replaced with 1...22 by the scripts.
export PHASE3_REF_PATTERN="${PHASE3_REF_PATTERN:-/path/to/1000G/chr{CHROM}.bcf}"
export PHASE3_SITES_TSV_PATTERN="${PHASE3_SITES_TSV_PATTERN:-/path/to/1000G/chr{CHROM}.sites.tsv.gz}"
export GENETIC_MAP_PATTERN="${GENETIC_MAP_PATTERN:-/path/to/maps/chr{CHROM}.gmap.gz}"
export AF_TABLE_PATTERN="${AF_TABLE_PATTERN:-/path/to/1000G_AF/v51.1_1240k_AF_ch{CHROM}.fixed.tsv}"

# Tools. Replace any name not available on PATH with an executable path.
export ADAPTERREMOVAL="${ADAPTERREMOVAL:-AdapterRemoval}"
export BWA="${BWA:-bwa}"
export SAMTOOLS="${SAMTOOLS:-samtools}"
export BCFTOOLS="${BCFTOOLS:-bcftools}"
export BGZIP="${BGZIP:-bgzip}"
export TABIX="${TABIX:-tabix}"
export ANGSD="${ANGSD:-angsd}"
export ANGSD_CONTAMINATION="${ANGSD_CONTAMINATION:-contamination}"
export TRIMBAM="${TRIMBAM:-bam}"
export MAPDAMAGE="${MAPDAMAGE:-mapDamage}"
export BAMREFINE="${BAMREFINE:-bamrefine}"
export GLIMPSE1_CHUNK="${GLIMPSE1_CHUNK:-GLIMPSE_chunk}"
export GLIMPSE1_PHASE="${GLIMPSE1_PHASE:-GLIMPSE_phase}"
export GLIMPSE1_LIGATE="${GLIMPSE1_LIGATE:-GLIMPSE_ligate}"
export GLIMPSE1_SAMPLE="${GLIMPSE1_SAMPLE:-GLIMPSE_sample}"
export GLIMPSE2_IMAGE="${GLIMPSE2_IMAGE:-glimpse:v2.0.0-27-g0919952_20221207}"
export DOCKER="${DOCKER:-docker}"
# Common host root mounted read-only at /data for GLIMPSE2 concordance.
export CONCORDANCE_INPUT_ROOT="${CONCORDANCE_INPUT_ROOT:-/path/to/concordance_inputs}"
export CONCORDANCE_FREQ_PATTERN="${CONCORDANCE_FREQ_PATTERN:-/path/to/concordance_inputs/frequency/chr{CHROM}.vcf.gz}"
export ANCIBD_RUN="${ANCIBD_RUN:-ancIBD-run}"
export ANCIBD_SUMMARY="${ANCIBD_SUMMARY:-ancIBD-summary}"
export PILEUPCALLER="${PILEUPCALLER:-pileupCaller}"
export READ2_PY="${READ2_PY:-/path/to/READv2/READ2.py}"
export PLINK="${PLINK:-plink}"
export PLINK2="${PLINK2:-plink2}"
export HAPLOGREP3="${HAPLOGREP3:-haplogrep3}"
export SMARTPCA="${SMARTPCA:-smartpca}"
export CONVERTF="${CONVERTF:-convertf}"

# schmutzi 1.5.7 executables
export BAM2PROF="${BAM2PROF:-bam2prof}"
export ENDOCALLER="${ENDOCALLER:-endoCaller}"
export MTCONT="${MTCONT:-mtCont}"
export CONTOUT2CONTEST="${CONTOUT2CONTEST:-contOut2ContEst.pl}"
