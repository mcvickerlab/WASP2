#!/bin/bash
#$ -N wasp2_rust_perf
#$ -t 1-150
#$ -tc 4
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# WASP2 Rust Performance Benchmark
# Based on aho's benchmark script

# Activate conda env
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Target subset
n_subset=$(($SGE_TASK_ID * 1000000))

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
log_file="${WASP2_DIR}/benchmarking/results/wasp2_rust_perf_${JOB_ID}.tsv"

genome_index="/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
input_bam="/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam"
input_vcf="/iblm/netapp/data1/aho/variants/NA12878.vcf.gz"
sample="NA12878"

# Total mapped reads in BAM
bam_reads=159133307

# Proportion to subsample
bam_prop=$(echo "scale=7; ${n_subset}/${bam_reads}" | bc)

# Random seed
seed=$RANDOM

# Create temp directory
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT

# Subset the BAM
prefix="bam_${n_subset}_${seed}"
subset_out="${dir}/${prefix}.bam"

samtools view -h --bam --subsample ${bam_prop} --subsample-seed ${seed} -o ${subset_out} ${input_bam}
samtools index ${subset_out}

# Process WASP2 with Rust
start=$(date +%s)

# IMPORTANT: Must use -m and PYTHONPATH for proper module import
export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"
python -m mapping make-reads ${subset_out} ${input_vcf} -s ${sample} --out ${dir} --paired --phased --temp_loc ${dir} --threads 8

intersect_end=$(date +%s)
intersect_runtime=$(($intersect_end-$start))

# Remap reads
r1_reads="${dir}/${prefix}_swapped_alleles_r1.fq"
r2_reads="${dir}/${prefix}_swapped_alleles_r2.fq"
remapped_bam="${dir}/${prefix}_remapped.bam"

bwa mem -t 16 -M ${genome_index} ${r1_reads} ${r2_reads} | samtools view -S -b -h -F 4 - > ${remapped_bam}
samtools sort -o ${remapped_bam} ${remapped_bam}
samtools index ${remapped_bam}

remap_end=$(date +%s)
remap_runtime=$(($remap_end-$intersect_end))

# Filter remapped reads
wasp_json="${dir}/${prefix}_wasp_data_files.json"
python -m mapping filter-remapped ${remapped_bam} --json ${wasp_json} --threads 8

end=$(date +%s)
filt_runtime=$(($end-$remap_end))
total_runtime=$(($end-$start))

# Append to log file
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${n_subset}" "${seed}" "$total_runtime" "$intersect_runtime" "$remap_runtime" "$filt_runtime" >> ${log_file}

echo "Task ${SGE_TASK_ID} complete: ${n_subset} reads in ${total_runtime}s"

# End-of-job summary
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
