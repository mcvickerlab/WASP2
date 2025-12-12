#!/bin/bash
# WASP1 (Python) benchmark on GM12878 ATAC-seq data - SNP only mode
# Full pipeline: BWA-MEM + WASP1 (find_intersecting_snps + filter_remapped_reads)
#
# For Figure 1B (timing) and Figure 1C (pre/post remap read counts)
#
#$ -N wasp1_gm12878_snp
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

# Timestamp for this run
TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "Benchmark timestamp: ${TIMESTAMP}"

# Paths
BENCHMARK_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878"
OUTPUT_DIR="${BENCHMARK_DIR}/results/wasp1_snp_${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"

# Data files
INPUT_BAM="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
REF_GENOME="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

# WASP1 source
WASP1_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/wasp1_source"
SNP2H5="${WASP1_DIR}/snp2h5/snp2h5"
FIND_INTERSECTING="${WASP1_DIR}/mapping/find_intersecting_snps.py"
FILTER_REMAPPED="${WASP1_DIR}/mapping/filter_remapped_reads.py"
CHROM_SIZES="${WASP1_DIR}/hg38_chromsizes.txt"

# Sample
SAMPLE="NA12878"
THREADS=8

# Conda
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Tools
BWA=$(which bwa)
SAMTOOLS=$(which samtools)
BEDTOOLS=$(which bedtools)
TIME="/usr/bin/time"

cd "${OUTPUT_DIR}"

# Log file
PROFILE_LOG="${OUTPUT_DIR}/profile.log"
echo "WASP1 SNP-only Pipeline Benchmark - GM12878 ATAC-seq" > ${PROFILE_LOG}
echo "======================================================" >> ${PROFILE_LOG}
echo "Timestamp: ${TIMESTAMP}" >> ${PROFILE_LOG}
echo "Date: $(date)" >> ${PROFILE_LOG}
echo "Sample: ${SAMPLE}" >> ${PROFILE_LOG}
echo "Threads: ${THREADS}" >> ${PROFILE_LOG}
echo "" >> ${PROFILE_LOG}

echo "========================================"
echo "WASP1 SNP Benchmark on GM12878 ATAC-seq"
echo "Timestamp: ${TIMESTAMP}"
echo "========================================"
echo "Start time: $(date)"
echo "Input BAM: ${INPUT_BAM}"
echo "VCF: ${VCF}"
echo "Threads: ${THREADS}"
echo ""

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: Copy input BAM to working directory (for output consistency)
# -----------------------------------------------------------------------------
echo "Copying input BAM to output directory..."
cp "${INPUT_BAM}" "${OUTPUT_DIR}/original.bam"
samtools index "${OUTPUT_DIR}/original.bam"
ORIGINAL_BAM="${OUTPUT_DIR}/original.bam"

# Count original reads
ORIGINAL_READS=$(samtools view -c ${ORIGINAL_BAM})
echo "Original reads: ${ORIGINAL_READS}"

# -----------------------------------------------------------------------------
# PRE-BENCHMARK: snp2h5 VCF -> HDF5 conversion
# -----------------------------------------------------------------------------
echo "========================================"
echo "PRE-BENCHMARK: snp2h5 VCF -> HDF5 conversion"
echo "========================================"

SNP2H5_DIR="${OUTPUT_DIR}/snp2h5_data"
mkdir -p "${SNP2H5_DIR}"

SNP2H5_PREP_START=$(date +%s.%N)

# Filter VCF to SNPs only and split by chromosome
echo "Filtering VCF to SNPs only and splitting by chromosome..."
for i in {1..22}; do
    echo "  Processing chr${i}..."
    bcftools view -r "chr${i}" -v snps ${VCF} -Oz -o ${SNP2H5_DIR}/chr${i}.vcf.gz 2>/dev/null || \
    bcftools view -r "${i}" -v snps ${VCF} -Oz -o ${SNP2H5_DIR}/chr${i}.vcf.gz 2>/dev/null || true
done

echo "Running snp2h5..."
${SNP2H5} --chrom ${CHROM_SIZES} \
          --format "vcf" \
          --haplotype ${SNP2H5_DIR}/haplotypes.h5 \
          --snp_index ${SNP2H5_DIR}/snp_index.h5 \
          --snp_tab ${SNP2H5_DIR}/snp_tab.h5 \
          ${SNP2H5_DIR}/*.vcf.gz 2>/dev/null || echo "snp2h5 completed with warnings"

SNP2H5_PREP_END=$(date +%s.%N)
SNP2H5_PREP_TIME=$(echo "${SNP2H5_PREP_END} - ${SNP2H5_PREP_START}" | bc)
echo "snp2h5 preprocessing time: ${SNP2H5_PREP_TIME}s (NOT counted in main benchmark)"

# Extract variant BED for Figure 1C analysis
echo "Creating variant BED file for read counting..."
bcftools query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' ${VCF} | \
    awk '$4 !~ /,/ && $5 !~ /,/ && length($4)==1 && length($5)==1' > ${OUTPUT_DIR}/variants.bed
VARIANT_COUNT=$(wc -l < ${OUTPUT_DIR}/variants.bed)
echo "SNP variants in BED: ${VARIANT_COUNT}"

# Count SNV-overlapping reads in original BAM (for Figure 1C pre-WASP)
echo "Counting SNV-overlapping reads in original BAM..."
SNV_READS_PRE=$(${BEDTOOLS} intersect -a ${ORIGINAL_BAM} -b ${OUTPUT_DIR}/variants.bed -u | \
    ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
echo "SNV-overlapping read pairs (pre-WASP): ${SNV_READS_PRE}"

# -----------------------------------------------------------------------------
# BENCHMARK STARTS HERE
# -----------------------------------------------------------------------------
echo "========================================"
echo "BENCHMARK TIMING STARTS NOW"
echo "========================================"

TOTAL_START=$(date +%s.%N)

# -----------------------------------------------------------------------------
# STEP 1: WASP1 find_intersecting_snps.py
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 1: WASP1 find_intersecting_snps.py"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 1: WASP1 find_intersecting_snps.py" >> ${PROFILE_LOG}

STEP1_START=$(date +%s.%N)

${TIME} -v python ${FIND_INTERSECTING} \
    --is_paired_end \
    --is_sorted \
    --output_dir ${OUTPUT_DIR} \
    --snp_index ${SNP2H5_DIR}/snp_index.h5 \
    --snp_tab ${SNP2H5_DIR}/snp_tab.h5 \
    --haplotype ${SNP2H5_DIR}/haplotypes.h5 \
    --samples ${SAMPLE} \
    original.bam 2>> ${PROFILE_LOG}

STEP1_END=$(date +%s.%N)
STEP1_TIME=$(echo "${STEP1_END} - ${STEP1_START}" | bc)
echo "STEP 1 completed in ${STEP1_TIME} seconds"

# Check for remap reads
REMAP_R1=$(ls original.remap.fq1.gz 2>/dev/null || echo "")
if [ -n "${REMAP_R1}" ]; then
    R1_READS=$(zcat original.remap.fq1.gz | wc -l)
    R1_PAIRS=$((R1_READS / 4))
    echo "Reads to remap: ${R1_PAIRS} pairs"
fi

# -----------------------------------------------------------------------------
# STEP 2: BWA-MEM remap flipped reads
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 2: BWA-MEM remap flipped reads"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 2: BWA-MEM Remap" >> ${PROFILE_LOG}

STEP2_START=$(date +%s.%N)

# Convert FASTQs to interleaved for BWA
zcat original.remap.fq1.gz > remap_r1.fq
zcat original.remap.fq2.gz > remap_r2.fq

${TIME} -v ${BWA} mem -t ${THREADS} ${REF_GENOME} remap_r1.fq remap_r2.fq 2>> ${PROFILE_LOG} | \
    ${SAMTOOLS} view -bS - | \
    ${SAMTOOLS} sort -@ ${THREADS} -o remapped.bam -

${SAMTOOLS} index remapped.bam
rm -f remap_r1.fq remap_r2.fq

STEP2_END=$(date +%s.%N)
STEP2_TIME=$(echo "${STEP2_END} - ${STEP2_START}" | bc)
echo "STEP 2 completed in ${STEP2_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 3: WASP1 filter_remapped_reads.py
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 3: WASP1 filter_remapped_reads.py"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 3: WASP1 filter_remapped_reads.py" >> ${PROFILE_LOG}

STEP3_START=$(date +%s.%N)

${TIME} -v python ${FILTER_REMAPPED} \
    original.to.remap.bam \
    remapped.bam \
    remap_filtered.bam 2>> ${PROFILE_LOG}

STEP3_END=$(date +%s.%N)
STEP3_TIME=$(echo "${STEP3_END} - ${STEP3_START}" | bc)
echo "STEP 3 completed in ${STEP3_TIME} seconds"

# -----------------------------------------------------------------------------
# STEP 4: Merge BAMs to create final keep BAM
# -----------------------------------------------------------------------------
echo "========================================"
echo "STEP 4: Merge BAMs (keep + filtered)"
echo "========================================"
echo "" >> ${PROFILE_LOG}
echo "STEP 4: Merge BAMs" >> ${PROFILE_LOG}

STEP4_START=$(date +%s.%N)

# The "remap_keep.bam" equivalent is the merged output
${SAMTOOLS} merge -@ ${THREADS} -o wasp1_merged.bam original.keep.bam remap_filtered.bam
${SAMTOOLS} sort -@ ${THREADS} -o remap_keep.bam wasp1_merged.bam
${SAMTOOLS} index remap_keep.bam
rm -f wasp1_merged.bam

STEP4_END=$(date +%s.%N)
STEP4_TIME=$(echo "${STEP4_END} - ${STEP4_START}" | bc)
echo "STEP 4 completed in ${STEP4_TIME} seconds"

# -----------------------------------------------------------------------------
# FIGURE 1C: SNV-specific read counting (post-WASP)
# -----------------------------------------------------------------------------
echo "========================================"
echo "FIGURE 1C: Counting SNV-overlapping reads"
echo "========================================"

# Count SNV-overlapping reads in filtered BAM (post-WASP)
echo "Counting SNV-overlapping reads in filtered BAM..."
SNV_READS_POST=$(${BEDTOOLS} intersect -a remap_keep.bam -b ${OUTPUT_DIR}/variants.bed -u | \
    ${SAMTOOLS} view | cut -f1 | sort -u | wc -l)
echo "SNV-overlapping read pairs (post-WASP): ${SNV_READS_POST}"

# Calculate SNV pass rate
if [ ${SNV_READS_PRE} -gt 0 ]; then
    SNV_PASS_RATE=$(echo "scale=2; ${SNV_READS_POST} * 100 / ${SNV_READS_PRE}" | bc)
else
    SNV_PASS_RATE=0
fi
echo "SNV read pass rate: ${SNV_PASS_RATE}%"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "${TOTAL_END} - ${TOTAL_START}" | bc)

# Count kept reads
KEPT_READS=$(samtools view -c remap_keep.bam)
PASS_RATE=$(echo "scale=2; ${KEPT_READS} * 100 / ${ORIGINAL_READS}" | bc)

# WASP-only time (excluding alignment)
WASP_ONLY=$(echo "${STEP1_TIME} + ${STEP3_TIME} + ${STEP4_TIME}" | bc)

echo "========================================"
echo "BENCHMARK RESULTS"
echo "========================================"
echo ""
echo "Individual step times:"
echo "  STEP 1 (WASP1 find_intersecting):  ${STEP1_TIME} s"
echo "  STEP 2 (BWA-MEM remap):            ${STEP2_TIME} s"
echo "  STEP 3 (WASP1 filter_remapped):    ${STEP3_TIME} s"
echo "  STEP 4 (Merge BAMs):               ${STEP4_TIME} s"
echo ""
echo "WASP-only time (steps 1+3+4):        ${WASP_ONLY} s"
echo "TOTAL time:                          ${TOTAL_TIME} s"
echo ""
echo "Read counts (total):"
echo "  Original reads:  ${ORIGINAL_READS}"
echo "  Kept reads:      ${KEPT_READS}"
echo "  Pass rate:       ${PASS_RATE}%"
echo ""
echo "SNV-overlapping reads (for Figure 1C):"
echo "  Pre-WASP:   ${SNV_READS_PRE} read pairs"
echo "  Post-WASP:  ${SNV_READS_POST} read pairs"
echo "  Pass rate:  ${SNV_PASS_RATE}%"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"

# Save results to JSON
cat > benchmark_results.json << EOF
{
    "timestamp": "${TIMESTAMP}",
    "pipeline": "wasp1_snp",
    "sample": "${SAMPLE}",
    "data_type": "ATAC-seq",
    "threads": ${THREADS},
    "step1_wasp1_find_intersecting_s": ${STEP1_TIME},
    "step2_bwa_remap_s": ${STEP2_TIME},
    "step3_wasp1_filter_remapped_s": ${STEP3_TIME},
    "step4_merge_bams_s": ${STEP4_TIME},
    "wasp_only_s": ${WASP_ONLY},
    "total_s": ${TOTAL_TIME},
    "snp2h5_preprocessing_s": ${SNP2H5_PREP_TIME},
    "original_reads": ${ORIGINAL_READS},
    "kept_reads": ${KEPT_READS},
    "pass_rate_percent": ${PASS_RATE},
    "snv_count": ${VARIANT_COUNT},
    "snv_reads_pre": ${SNV_READS_PRE},
    "snv_reads_post": ${SNV_READS_POST},
    "snv_pass_rate_percent": ${SNV_PASS_RATE}
}
EOF

echo ""
echo "Results saved to: ${OUTPUT_DIR}/benchmark_results.json"
echo "Key output files for Figure 1C:"
echo "  - original.bam (pre-remap)"
echo "  - remap_keep.bam (post-remap)"
echo "  - variants.bed (SNP positions)"
