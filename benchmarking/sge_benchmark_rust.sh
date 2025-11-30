#!/bin/bash
#$ -N wasp2_rust_perf
#$ -t 1-150
#$ -tc 4
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -e /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# WASP2 Rust Performance Benchmark
# Tests scaling from 1M to 150M reads with 8 threads

source ~/.bashrc
conda activate wasp2

N_READS=$((SGE_TASK_ID * 1000000))
SEED=$RANDOM

WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
DATA_DIR="/iblm/netapp/home/aho/projects/wasp/testing/performance/data"
OUTPUT_DIR="${WASP2_DIR}/benchmarking/results"

BAM_FILE="${DATA_DIR}/CVPC_ATAC_chr10_30x.bam"
VCF_FILE="${DATA_DIR}/CVPC_ATAC_chr10.vcf.gz"
REF_FASTA="/iblm/netapp/data1/external/GRCh38/GRCh38.fa"

TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Subsample BAM
SUBSAMPLE_BAM="${TEMP_DIR}/subsample.bam"
FRAC=$(echo "scale=6; ${N_READS} / 30000000" | bc)
if (( $(echo "$FRAC > 1" | bc -l) )); then
    FRAC=1.0
fi
samtools view -@ 8 -s ${SEED}.${FRAC} -b ${BAM_FILE} > ${SUBSAMPLE_BAM}
samtools index ${SUBSAMPLE_BAM}

# Run WASP2 with Rust and time each stage
START_TOTAL=$(date +%s)

# Intersect stage
START_INTERSECT=$(date +%s)
python -m mapping intersect \
    ${SUBSAMPLE_BAM} \
    ${VCF_FILE} \
    --output-dir ${TEMP_DIR} \
    --threads 8
END_INTERSECT=$(date +%s)

# Remap stage
START_REMAP=$(date +%s)
python -m mapping remap \
    ${TEMP_DIR}/to_remap.bam \
    ${REF_FASTA} \
    --output-dir ${TEMP_DIR} \
    --threads 8
END_REMAP=$(date +%s)

# Filter stage
START_FILTER=$(date +%s)
python -m mapping filter \
    ${SUBSAMPLE_BAM} \
    ${TEMP_DIR}/remapped.bam \
    --output-dir ${TEMP_DIR} \
    --threads 8
END_FILTER=$(date +%s)

END_TOTAL=$(date +%s)

# Calculate times
TIME_TOTAL=$((END_TOTAL - START_TOTAL))
TIME_INTERSECT=$((END_INTERSECT - START_INTERSECT))
TIME_REMAP=$((END_REMAP - START_REMAP))
TIME_FILTER=$((END_FILTER - START_FILTER))

# Write results
RESULTS_FILE="${OUTPUT_DIR}/wasp2_rust_perf_${JOB_ID}.tsv"
echo -e "${N_READS}\t${SEED}\t${TIME_TOTAL}\t${TIME_INTERSECT}\t${TIME_REMAP}\t${TIME_FILTER}" >> ${RESULTS_FILE}

echo "Task ${SGE_TASK_ID} complete: ${N_READS} reads in ${TIME_TOTAL}s"
