#!/bin/bash
set -e

# Quick test to verify --include-indels flag works
# Uses only chr22 from VCF (smallest chromosome)

cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

export PATH="/iblm/netapp/home/jjaureguy/mambaforge/envs/WASP2_new/bin:$PATH"
export PYTHONPATH="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/src:$PYTHONPATH"

BAM="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
GTF="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"

TESTDIR="gm12878_benchmark/test_quick"
mkdir -p ${TESTDIR}

echo "=========================================="
echo "Quick Test: --include-indels flag"
echo "=========================================="

# Create tiny GTF with just 3 genes for testing
echo "Creating test GTF with 3 imprinted genes (H19, XIST, IGF2)..."
grep -E "H19|XIST|IGF2" ${GTF} | head -100 > ${TESTDIR}/test_3genes.gtf

echo ""
echo "Test 1: SNPs only (no --include-indels)"
echo "----------------------------------------"
/iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.counting count-variants \
    ${BAM} ${VCF} \
    -s NA12878 \
    -r ${TESTDIR}/test_3genes.gtf \
    -o ${TESTDIR}/counts_SNPs_ONLY.tsv \
    --gene_feature transcript \
    --gene_attribute transcript_id \
    --gene_parent gene_name \
    2>&1 | head -20

SNP_COUNT=$(wc -l < ${TESTDIR}/counts_SNPs_ONLY.tsv)
echo "✅ SNPs only: ${SNP_COUNT} variants"

echo ""
echo "Test 2: SNPs + indels (WITH --include-indels)"
echo "----------------------------------------------"
/iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.counting count-variants \
    ${BAM} ${VCF} \
    -s NA12878 \
    -r ${TESTDIR}/test_3genes.gtf \
    -o ${TESTDIR}/counts_WITH_INDELS.tsv \
    --gene_feature transcript \
    --gene_attribute transcript_id \
    --gene_parent gene_name \
    --include-indels \
    2>&1 | head -20

INDEL_COUNT=$(wc -l < ${TESTDIR}/counts_WITH_INDELS.tsv)
echo "✅ SNPs+indels: ${INDEL_COUNT} variants"

echo ""
echo "=========================================="
echo "RESULTS"
echo "=========================================="
echo "SNPs only:     ${SNP_COUNT} variants"
echo "SNPs + indels: ${INDEL_COUNT} variants"

if [ ${INDEL_COUNT} -gt ${SNP_COUNT} ]; then
    DIFF=$((INDEL_COUNT - SNP_COUNT))
    echo "New indels:    ${DIFF} variants"
    echo ""
    echo "✅ SUCCESS: --include-indels flag is working!"
    echo "   Indel count increased as expected."
else
    echo ""
    echo "⚠️  WARNING: No increase in variant count!"
    echo "   The --include-indels flag may not be working."
fi

echo ""
echo "Quick check of variant types:"
echo "---"
head -5 ${TESTDIR}/counts_WITH_INDELS.tsv
echo "..."
echo ""
echo "Test files saved in: ${TESTDIR}/"
