#!/bin/bash
set -e

# Go to project root directory
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

# Add bioinformatics tools to PATH (bedtools, bcftools, etc.)
export PATH="/iblm/netapp/home/jjaureguy/mambaforge/envs/WASP2_new/bin:$PATH"

# Set PYTHONPATH to include src directory
export PYTHONPATH="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/src:$PYTHONPATH"

# Use base Python 3.10 which has all dependencies installed

BAM="/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
VCF="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
GTF="/iblm/netapp/data3/aho/project_data/wasp2/imprinted_rna/data/geneimprint.gtf"

OUTDIR="gm12878_benchmark/results"
LOGDIR="gm12878_benchmark/logs"

echo "Step 1: Counting variants (SNPs + indels)..."
/iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.counting count-variants \
    ${BAM} ${VCF} \
    -s NA12878 \
    -r ${GTF} \
    -o ${OUTDIR}/counts_WITH_INDELS.tsv \
    --gene_feature transcript \
    --gene_attribute transcript_id \
    --gene_parent gene_name \
    --include-indels \
    2>&1 | tee ${LOGDIR}/counting.log

echo "Step 2: Finding allelic imbalance..."
/iblm/netapp/home/jjaureguy/mambaforge/bin/python -m src.analysis find-imbalance \
    --phased \
    --out ${OUTDIR}/ai_results_WITH_INDELS.tsv \
    --group gene_name \
    ${OUTDIR}/counts_WITH_INDELS.tsv \
    2>&1 | tee ${LOGDIR}/analysis.log

echo "Done! Results in ${OUTDIR}/"
