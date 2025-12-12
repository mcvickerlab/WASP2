#!/bin/bash
# Prepare chr22 subset for WASP2-Rust vs Python DEV comparison
#
#$ -N prepare_chr22
#$ -V
#$ -pe iblm 4
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Configuration
INPUT_BAM="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
INPUT_VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
REFERENCE="/iblm/netapp/data1/external/GRCh38_no_alt/genome.fa"

# Output directory
WORKDIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison"
mkdir -p "${WORKDIR}"/{input,rust,python,analysis}
mkdir -p "${WORKDIR}"/rust/{step1_vcf_to_bed,step2_filter_bam,step3_make_reads,step4_remap,step5_filter_remapped}
mkdir -p "${WORKDIR}"/python/{step1_vcf_to_bed,step2_filter_bam,step3_make_reads,step4_remap,step5_filter_remapped}

cd "${WORKDIR}"

echo "========================================"
echo "Preparing chr22 Subset"
echo "Timestamp: $(date)"
echo "========================================"

# 1. Extract chr22 from BAM
echo ""
echo "Step 1: Extracting chr22 from BAM..."
if [[ ! -f input/chr22_input.bam ]]; then
    samtools view -@ 4 -b "${INPUT_BAM}" chr22 > input/chr22_input.bam
    samtools index -@ 4 input/chr22_input.bam
    echo "  chr22 BAM reads: $(samtools view -c input/chr22_input.bam)"
else
    echo "  chr22 BAM already exists: $(samtools view -c input/chr22_input.bam) reads"
fi

# 2. Extract chr22 from VCF
echo ""
echo "Step 2: Extracting chr22 from VCF..."
if [[ ! -f input/chr22_variants.vcf.gz ]]; then
    bcftools view -r chr22 "${INPUT_VCF}" -Oz -o input/chr22_variants.vcf.gz
    tabix -p vcf input/chr22_variants.vcf.gz
    echo "  chr22 variants: $(bcftools view -H input/chr22_variants.vcf.gz | wc -l)"
else
    echo "  chr22 VCF already exists: $(bcftools view -H input/chr22_variants.vcf.gz | wc -l) variants"
fi

# 3. Create symlink to reference
echo ""
echo "Step 3: Setting up reference..."
ln -sf "${REFERENCE}" input/reference.fa
ln -sf "${REFERENCE}.fai" input/reference.fa.fai

# Find BWA index files
BWA_PREFIX="${REFERENCE%.*}"
for ext in amb ann bwt pac sa; do
    if [[ -f "${BWA_PREFIX}.${ext}" ]]; then
        ln -sf "${BWA_PREFIX}.${ext}" input/reference.${ext}
    elif [[ -f "${REFERENCE}.${ext}" ]]; then
        ln -sf "${REFERENCE}.${ext}" input/reference.fa.${ext}
    fi
done

# Summary
echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo "Working directory: ${WORKDIR}"
echo "chr22 BAM: input/chr22_input.bam"
echo "chr22 VCF: input/chr22_variants.vcf.gz"
echo "Reference: input/reference.fa"
echo ""
echo "Ready for pipeline comparison scripts."
echo "Completed: $(date)"
