#!/bin/bash
#$ -N wasp2_rust_indels
#$ -t 1-150
#$ -tc 16
#$ -V
#$ -pe iblm 8
#$ -l h_vmem=8G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# WASP2 Rust Performance Benchmark WITH INDELS
# Tests read counts from 1M to 150M in 1M increments (150 tasks) - matches SNP benchmark

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Threading controls
# - THREADS defaults to NSLOTS (SGE), else 8
# - COMPRESSION_THREADS is per FASTQ file (R1 and R2); defaults to 1 to avoid oversubscription
THREADS="${THREADS:-${NSLOTS:-8}}"
COMPRESSION_THREADS="${COMPRESSION_THREADS:-1}"

# Read counts: 1M, 2M, 3M, ..., 150M (150 tasks)
n_subset=$(($SGE_TASK_ID * 1000000))

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
log_file="${WASP2_DIR}/benchmarking/results/wasp2_rust_indels_scale_${JOB_ID}.tsv"

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

# Process WASP2 with Rust (WITH INDELS)
start=$(date +%s)

export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"

# Use unified pipeline with indels enabled
python -c "
from mapping.run_mapping import run_make_remap_reads_unified

	stats = run_make_remap_reads_unified(
	    bam_file='${subset_out}',
	    variant_file='${input_vcf}',
	    samples='${sample}',
	    out_dir='${dir}',
	    include_indels=True,
	    max_indel_len=10,
	    threads=${THREADS},
	    compression_threads=${COMPRESSION_THREADS},
	    use_parallel=True
	)
	print(f'Pairs processed: {stats[\"pairs_processed\"]:,}')
	print(f'Haplotypes written: {stats[\"haplotypes_written\"]:,}')
"

intersect_end=$(date +%s)
intersect_runtime=$(($intersect_end-$start))

# Find FASTQ files
r1_reads=$(ls ${dir}/*_remap_r1.fq.gz 2>/dev/null | head -1)
r2_reads=$(ls ${dir}/*_remap_r2.fq.gz 2>/dev/null | head -1)

# Remap reads
zcat ${r1_reads} > ${dir}/remap_r1.fq
zcat ${r2_reads} > ${dir}/remap_r2.fq
remapped_bam="${dir}/${prefix}_remapped.bam"

bwa mem -t 16 -M ${genome_index} ${dir}/remap_r1.fq ${dir}/remap_r2.fq | samtools view -S -b -h -F 4 - > ${remapped_bam}
samtools sort -o ${remapped_bam} ${remapped_bam}
samtools index ${remapped_bam}

remap_end=$(date +%s)
remap_runtime=$(($remap_end-$intersect_end))

# Filter remapped reads
python -c "
from wasp2_rust import filter_bam_wasp

	kept, removed_moved, removed_missing = filter_bam_wasp(
	    '${subset_out}',
	    '${remapped_bam}',
	    '${dir}/remap_keep.bam',
	    threads=${THREADS}
	)
	print(f'Kept: {kept}, Removed: {removed_moved}, Missing: {removed_missing}')
	"

end=$(date +%s)
filt_runtime=$(($end-$remap_end))
total_runtime=$(($end-$start))

# Append to log file
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${n_subset}" "${seed}" "$total_runtime" "$intersect_runtime" "$remap_runtime" "$filt_runtime" >> ${log_file}

echo "Task ${SGE_TASK_ID} complete: ${n_subset} reads (with INDELS) in ${total_runtime}s"

# End-of-job summary
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
