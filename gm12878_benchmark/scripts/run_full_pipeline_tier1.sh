#!/bin/bash
set -e

# Tier 1: Quick Full Pipeline Test with Indels
# Goal: Verify all 4 steps work end-to-end
# Scope: 3 genes (H19, IGF2, SNRPN)
# Expected runtime: <5 minutes

cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

export PATH="/iblm/netapp/home/jjaureguy/mambaforge/envs/WASP2_new/bin:$PATH"
export PYTHONPATH="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/src:$PYTHONPATH"

# Enable Rust acceleration
export WASP2_DISABLE_RUST=0

TESTDIR="gm12878_benchmark/test_full_pipeline"
mkdir -p ${TESTDIR}/logs

THREADS="${THREADS:-8}"
BWA_THREADS="${BWA_THREADS:-${THREADS}}"
SAMTOOLS_SORT_THREADS="${SAMTOOLS_SORT_THREADS:-${THREADS}}"

BAM="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
GTF="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"
REF="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

echo "=========================================="
echo "TIER 1: Full WASP2 Pipeline Test (3 genes)"
echo "=========================================="
echo "Testing indel support across entire pipeline"
echo ""

# Create 3-gene GTF
echo "Creating test GTF with 3 imprinted genes (H19, IGF2, SNRPN)..."
grep -E "H19|IGF2|SNRPN" ${GTF} | head -100 > ${TESTDIR}/test_3genes.gtf
echo "✅ Test GTF created"
echo ""

# ==========================================
# STEP 1: Generate alternative reads
# ==========================================
echo "Step 1: Generating alternative reads with indels..."
echo "----------------------------------------"
/usr/bin/time -v /iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.mapping make-reads \
    ${BAM} ${VCF} \
    -s NA12878 \
    --indels \
    --max-indel-len 10 \
    --insert-qual 30 \
    --max-seqs 64 \
    -o ${TESTDIR}/remap_data \
    2>&1 | tee ${TESTDIR}/logs/step1_make_reads.log

if [ -f "${TESTDIR}/remap_data/to_remap.fq1.gz" ]; then
    READ_COUNT=$(zcat ${TESTDIR}/remap_data/to_remap.fq1.gz | wc -l)
    READ_COUNT=$((READ_COUNT / 4))
    echo "✅ Step 1 complete: ${READ_COUNT} read pairs generated"
else
    echo "❌ Step 1 failed: No FASTQ files generated"
    exit 1
fi
echo ""

# ==========================================
# STEP 2: Remap with BWA
# ==========================================
echo "Step 2: Remapping with BWA (${BWA_THREADS} threads)..."
echo "----------------------------------------"
/usr/bin/time -v bwa mem -t ${BWA_THREADS} \
    ${REF} \
    ${TESTDIR}/remap_data/to_remap.fq1.gz \
    ${TESTDIR}/remap_data/to_remap.fq2.gz \
    2> ${TESTDIR}/logs/step2_bwa_stderr.log \
    | samtools sort -@ ${SAMTOOLS_SORT_THREADS} -o ${TESTDIR}/remapped.bam - \
    2>&1 | tee ${TESTDIR}/logs/step2_bwa_remap.log

samtools index ${TESTDIR}/remapped.bam

if [ -f "${TESTDIR}/remapped.bam" ]; then
    REMAP_COUNT=$(samtools view -c ${TESTDIR}/remapped.bam)
    echo "✅ Step 2 complete: ${REMAP_COUNT} reads remapped"
else
    echo "❌ Step 2 failed: Remapped BAM not created"
    exit 1
fi
echo ""

# ==========================================
# STEP 3: Filter for concordant mappings
# ==========================================
echo "Step 3: Filtering remapped reads (with indel tolerance)..."
echo "----------------------------------------"
/usr/bin/time -v /iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.mapping filter-remapped \
    ${TESTDIR}/remapped.bam \
    --json ${TESTDIR}/remap_data/wasp_data_files.json \
    --same-locus-slop 2 \
    --use-rust \
    --threads 4 \
    -o ${TESTDIR}/wasp_filtered.bam \
    2>&1 | tee ${TESTDIR}/logs/step3_filter.log

if [ -f "${TESTDIR}/wasp_filtered.bam" ]; then
    FILT_COUNT=$(samtools view -c ${TESTDIR}/wasp_filtered.bam)
    echo "✅ Step 3 complete: ${FILT_COUNT} reads kept after WASP filtering"

    if [ ${REMAP_COUNT} -gt 0 ]; then
        CONCORDANCE=$(awk "BEGIN {printf \"%.1f\", 100 * ${FILT_COUNT} / ${REMAP_COUNT}}")
        echo "   Concordance rate: ${CONCORDANCE}%"
    fi
else
    echo "❌ Step 3 failed: Filtered BAM not created"
    exit 1
fi
echo ""

# ==========================================
# STEP 4: Count alleles
# ==========================================
echo "Step 4: Counting variants from WASP-filtered BAM..."
echo "----------------------------------------"
/usr/bin/time -v /iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.counting count-variants \
    ${TESTDIR}/wasp_filtered.bam ${VCF} \
    -s NA12878 \
    -r ${TESTDIR}/test_3genes.gtf \
    -o ${TESTDIR}/counts_WITH_INDELS_WASP_FILTERED.tsv \
    --gene_feature transcript \
    --gene_attribute transcript_id \
    --gene_parent gene_name \
    --include-indels \
    2>&1 | tee ${TESTDIR}/logs/step4_counting.log

if [ -f "${TESTDIR}/counts_WITH_INDELS_WASP_FILTERED.tsv" ]; then
    VAR_COUNT=$(wc -l < ${TESTDIR}/counts_WITH_INDELS_WASP_FILTERED.tsv)
    echo "✅ Step 4 complete: ${VAR_COUNT} variants counted"
else
    echo "❌ Step 4 failed: Counts file not created"
    exit 1
fi
echo ""

# ==========================================
# Summary
# ==========================================
echo "=========================================="
echo "TIER 1 TEST COMPLETE!"
echo "=========================================="
echo ""
echo "Pipeline Steps:"
echo "  1. Make reads:      ${READ_COUNT} read pairs generated"
echo "  2. BWA remap:       ${REMAP_COUNT} reads remapped"
echo "  3. WASP filter:     ${FILT_COUNT} reads kept (${CONCORDANCE}% concordance)"
echo "  4. Count variants:  ${VAR_COUNT} variants"
echo ""
echo "Check indel examples:"
head -20 ${TESTDIR}/counts_WITH_INDELS_WASP_FILTERED.tsv
echo "..."
echo ""
echo "Memory/runtime metrics saved in: ${TESTDIR}/logs/"
echo ""
echo "Extract metrics with:"
echo "  grep 'Maximum resident set size' ${TESTDIR}/logs/*.log"
echo "  grep 'Elapsed (wall clock)' ${TESTDIR}/logs/*.log"
echo ""
echo "Next: Review results, then run Tier 2 (chr22 test)"
