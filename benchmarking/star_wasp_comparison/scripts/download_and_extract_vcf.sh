#!/bin/bash
# Download 1000 Genomes VCFs and extract HG00731 heterozygous variants

set -euo pipefail

DATA_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/star_wasp_comparison/data"
VCF_DIR="${DATA_DIR}/1000g_vcf"
OUTPUT="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"
RENAME_CHRS="${DATA_DIR}/chr_rename.txt"

mkdir -p ${VCF_DIR}
cd ${VCF_DIR}

FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"

# Download all chromosomes (skip if already exists)
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    VCF="ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    TBI="${VCF}.tbi"

    if [ ! -f "${VCF}" ]; then
        echo "Downloading ${VCF}..."
        wget -c "${FTP_BASE}/${VCF}"
        wget -c "${FTP_BASE}/${TBI}"
    else
        echo "${VCF} already exists, skipping..."
    fi
done

echo "All VCFs downloaded."

# Extract HG00731 heterozygous variants
echo "Extracting HG00731 heterozygous variants..."

# Create list of VCF files
VCF_LIST=""
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    VCF_LIST="${VCF_LIST} ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
done

# Extract HG00731 and filter for heterozygous only, then rename contigs to `chr*`
# to match our STAR index/BAM reference names.
bcftools concat ${VCF_LIST} | \
    bcftools view -s HG00731 | \
    bcftools view -g het -Ou | \
    bcftools annotate --rename-chrs "${RENAME_CHRS}" -Oz -o ${OUTPUT}

bcftools index -t ${OUTPUT}

echo "Done! Output: ${OUTPUT}"
echo "Variant count:"
bcftools view -H ${OUTPUT} | wc -l
