#!/bin/bash
#$ -N wasp2_rust_indels
#$ -t 1-150
#$ -tc 10
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# =============================================================================
# WASP2 RUST INDEL SCALING BENCHMARK - 1M to 150M reads
# Same as SNP benchmark but with --indels flag enabled
# =============================================================================
# Task ID determines read count: task 1 = 1M, task 150 = 150M
#
# Threading (same as aho):
#   - SGE: 8 cores (-pe iblm 8)
#   - make-reads: -t 8
#   - BWA mem: -t 16
#   - filter-remapped: -t 8
#
# INDEL-specific Rust optimizations:
#   - VCF→BED with indels: Rust/noodles (include_indels=True)
#   - BAM-BED intersection: Rust/coitrees (variable-length variants)
#   - Remapping: Rust with CIGAR-aware position mapping
#   - Filtering: Rust with slop tolerance for micro-homology
# =============================================================================

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
log_file="${WASP2_DIR}/benchmarking/results/wasp2_rust_indels_${JOB_ID}.tsv"

# Same inputs as aho/SNP benchmark
genome_index="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
input_bam="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
input_vcf="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
sample="NA12878"

# Read count scales with task ID (same as aho)
n_subset=$(($SGE_TASK_ID * 1000000))
bam_reads=159133307
bam_prop=$(echo "scale=7; ${n_subset}/${bam_reads}" | bc)
seed=$RANDOM

dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT

echo "=============================================="
echo "WASP2 RUST INDEL SCALING - Task ${SGE_TASK_ID}"
echo "Reads: ${n_subset} (${SGE_TASK_ID}M)"
echo "INDELS: ENABLED"
echo "Date: $(date)"
echo "=============================================="

# Subset BAM
prefix="bam_${n_subset}_${seed}"
subset_bam="${dir}/${prefix}.bam"
samtools view -h --bam --subsample ${bam_prop} --subsample-seed ${seed} -o ${subset_bam} ${input_bam}
samtools index ${subset_bam}

# WASP2 make-reads with INDELS enabled (Rust accelerated + CIGAR-aware)
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"
start=$(date +%s)

python -m mapping make-reads ${subset_bam} ${input_vcf} -t 8 -s ${sample} --out ${dir} --paired --phased --indels --temp_loc ${dir}

intersect_end=$(date +%s)
intersect_runtime=$(($intersect_end-$start))

# BWA remap
r1_reads="${dir}/${prefix}_swapped_alleles_r1.fq"
r2_reads="${dir}/${prefix}_swapped_alleles_r2.fq"
remapped_bam="${dir}/${prefix}_remapped.bam"

bwa mem -t 16 -M ${genome_index} ${r1_reads} ${r2_reads} | samtools view -S -b -h -F 4 - > ${remapped_bam}
samtools sort -o ${remapped_bam} ${remapped_bam}
samtools index ${remapped_bam}

remap_end=$(date +%s)
remap_runtime=$(($remap_end-$intersect_end))

# WASP2 filter with slop for indel micro-homology (Rust accelerated)
wasp_json="${dir}/${prefix}_wasp_data_files.json"
python -m mapping filter-remapped ${remapped_bam} --json ${wasp_json} -t 8

end=$(date +%s)
filt_runtime=$(($end-$remap_end))
total_runtime=$(($end-$start))

echo ""
echo "=============================================="
echo "RESULTS - ${SGE_TASK_ID}M reads (INDELS)"
echo "=============================================="
echo "  Make-reads:  ${intersect_runtime}s"
echo "  BWA remap:   ${remap_runtime}s"
echo "  Filter:      ${filt_runtime}s"
echo "  Total:       ${total_runtime}s"
echo "=============================================="

# Same format as SNP benchmark: n_reads, seed, total, intersect, remap, filter
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${n_subset}" "${seed}" "$total_runtime" "$intersect_runtime" "$remap_runtime" "$filt_runtime" >> ${log_file}

echo "Task ${SGE_TASK_ID} complete (INDELS)"
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
