#!/bin/bash
# Run WASP2-Python DEV pipeline step-by-step on chr22 subset
#
#$ -N run_python_chr22
#$ -V
#$ -pe iblm 4
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

WORKDIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison"
PYTHON_DIR="${WORKDIR}/python"
INPUT_BAM="${WORKDIR}/input/chr22_input.bam"
INPUT_VCF="${WORKDIR}/input/chr22_variants.vcf.gz"
REFERENCE="${WORKDIR}/input/reference.fa"

cd "${WORKDIR}"

echo "========================================"
echo "WASP2-Python DEV Pipeline (chr22 subset)"
echo "Timestamp: $(date)"
echo "========================================"

# Python script for Python DEV module calls
python3 << 'PYTHON_SCRIPT'
import sys
import time
import json
import subprocess

# Add Python DEV source to path
sys.path.insert(0, '/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-python-dev/src/mapping')

from intersect_variant_data import vcf_to_bed, process_bam, intersect_reads
from make_remap_reads import write_remap_bam
from filter_remap_reads import filt_remapped_reads

WORKDIR = "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison"
PYTHON_DIR = f"{WORKDIR}/python"
INPUT_BAM = f"{WORKDIR}/input/chr22_input.bam"
INPUT_VCF = f"{WORKDIR}/input/chr22_variants.vcf.gz"
REFERENCE = f"{WORKDIR}/input/reference.fa"

results = {"pipeline": "wasp2_python_dev", "steps": {}}

print("\n" + "="*50)
print("STEP 1: VCF to BED")
print("="*50)
start = time.time()
bed_file = f"{PYTHON_DIR}/step1_vcf_to_bed/variants.bed"
vcf_to_bed(INPUT_VCF, bed_file, samples=["NA12878"])
elapsed = time.time() - start
results["steps"]["step1_vcf_to_bed"] = {"time_s": elapsed}

# Count variants
with open(bed_file) as f:
    variant_count = sum(1 for _ in f)
results["steps"]["step1_vcf_to_bed"]["variant_count"] = variant_count
print(f"  Variants in BED: {variant_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 2: Filter BAM by Variants (process_bam)")
print("="*50)
start = time.time()
to_remap_bam = f"{PYTHON_DIR}/step2_filter_bam/to_remap.bam"
keep_bam = f"{PYTHON_DIR}/step2_filter_bam/keep.bam"
remap_reads_file = f"{PYTHON_DIR}/step2_filter_bam/remap_reads.txt"
process_bam(INPUT_BAM, bed_file, to_remap_bam, remap_reads_file, keep_bam, threads=4, is_paired=True)
elapsed = time.time() - start
results["steps"]["step2_filter_bam"] = {"time_s": elapsed}

# Count reads
to_remap_count = int(subprocess.check_output(f"samtools view -c {to_remap_bam}", shell=True).decode().strip())
keep_count = int(subprocess.check_output(f"samtools view -c {keep_bam}", shell=True).decode().strip())
results["steps"]["step2_filter_bam"]["to_remap_reads"] = to_remap_count
results["steps"]["step2_filter_bam"]["keep_reads"] = keep_count

# Copy read names for comparison (already saved as remap_reads_file)
subprocess.run(f"cp {remap_reads_file} {PYTHON_DIR}/step2_filter_bam/read_names.txt", shell=True)
unique_names = int(subprocess.check_output(f"wc -l < {remap_reads_file}", shell=True).decode().strip())
results["steps"]["step2_filter_bam"]["unique_read_names"] = unique_names

print(f"  To remap reads: {to_remap_count}")
print(f"  Keep reads: {keep_count}")
print(f"  Unique read names to remap: {unique_names}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 2b: Intersect reads with variants")
print("="*50)
start = time.time()
intersect_file = f"{PYTHON_DIR}/step2_filter_bam/intersect.bed"
intersect_reads(to_remap_bam, bed_file, intersect_file)
elapsed = time.time() - start
results["steps"]["step2b_intersect"] = {"time_s": elapsed}

intersect_count = int(subprocess.check_output(f"wc -l < {intersect_file}", shell=True).decode().strip())
results["steps"]["step2b_intersect"]["intersect_count"] = intersect_count
print(f"  Intersection records: {intersect_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 3: Make Remap Reads (write_remap_bam)")
print("="*50)
start = time.time()
remap_fq1 = f"{PYTHON_DIR}/step3_make_reads/to_remap_R1.fq.gz"
remap_fq2 = f"{PYTHON_DIR}/step3_make_reads/to_remap_R2.fq.gz"
write_remap_bam(to_remap_bam, intersect_file, remap_fq1, remap_fq2, samples=["NA12878"], threads=4)
elapsed = time.time() - start
results["steps"]["step3_make_reads"] = {"time_s": elapsed}

# Count FASTQ entries (R1 count = pair count)
r1_count = int(subprocess.check_output(f"zcat {remap_fq1} | grep -c '^@'", shell=True).decode().strip())
r2_count = int(subprocess.check_output(f"zcat {remap_fq2} | grep -c '^@'", shell=True).decode().strip())
results["steps"]["step3_make_reads"]["r1_entries"] = r1_count
results["steps"]["step3_make_reads"]["r2_entries"] = r2_count
results["steps"]["step3_make_reads"]["fastq_entries"] = r1_count + r2_count

# Extract WASP names for analysis
subprocess.run(f"zcat {remap_fq1} | grep '^@' | sed 's/^@//' | cut -d' ' -f1 > {PYTHON_DIR}/step3_make_reads/wasp_names.txt", shell=True)

print(f"  R1 FASTQ entries: {r1_count}")
print(f"  R2 FASTQ entries: {r2_count}")
print(f"  Total FASTQ entries: {r1_count + r2_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 4: BWA-MEM Remap")
print("="*50)
start = time.time()
remapped_bam = f"{PYTHON_DIR}/step4_remap/remapped.bam"
# Python DEV uses separate R1/R2 files
subprocess.run(f"bwa mem -t 4 {REFERENCE} {remap_fq1} {remap_fq2} 2>/dev/null | samtools sort -@ 4 -o {remapped_bam}", shell=True, check=True)
subprocess.run(f"samtools index {remapped_bam}", shell=True, check=True)
elapsed = time.time() - start
results["steps"]["step4_remap"] = {"time_s": elapsed}

remapped_count = int(subprocess.check_output(f"samtools view -c {remapped_bam}", shell=True).decode().strip())
results["steps"]["step4_remap"]["remapped_reads"] = remapped_count
print(f"  Remapped reads: {remapped_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 5: Filter Remapped (filt_remapped_reads)")
print("="*50)
start = time.time()
keep_remapped = f"{PYTHON_DIR}/step5_filter_remapped/keep.bam"
keep_names_file = f"{PYTHON_DIR}/step5_filter_remapped/kept_names.txt"
filt_remapped_reads(to_remap_bam, remapped_bam, keep_remapped, keep_read_file=keep_names_file, threads=4)
elapsed = time.time() - start
results["steps"]["step5_filter_remapped"] = {"time_s": elapsed}

kept_count = int(subprocess.check_output(f"samtools view -c {keep_remapped}", shell=True).decode().strip())
results["steps"]["step5_filter_remapped"]["kept_reads"] = kept_count

# Calculate removed
removed_count = remapped_count - kept_count
results["steps"]["step5_filter_remapped"]["removed_reads"] = removed_count

# Copy kept names
subprocess.run(f"cp {keep_names_file} {PYTHON_DIR}/step5_filter_remapped/kept_read_names.txt", shell=True)

print(f"  Kept reads: {kept_count}")
print(f"  Removed reads: {removed_count}")
print(f"  Time: {elapsed:.2f}s")

# Save results
with open(f"{PYTHON_DIR}/pipeline_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*50)
print("PYTHON DEV PIPELINE COMPLETE")
print("="*50)
print(f"Results saved to: {PYTHON_DIR}/pipeline_results.json")

PYTHON_SCRIPT

echo ""
echo "Completed: $(date)"
