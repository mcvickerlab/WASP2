#!/bin/bash
# Run WASP2-Rust pipeline step-by-step on chr22 subset
#
#$ -N run_rust_chr22
#$ -V
#$ -pe iblm 4
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/logs/
#$ -cwd

set -e

WORKDIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison"
RUST_DIR="${WORKDIR}/rust"
INPUT_BAM="${WORKDIR}/input/chr22_input.bam"
INPUT_VCF="${WORKDIR}/input/chr22_variants.vcf.gz"
REFERENCE="${WORKDIR}/input/reference.fa"

cd "${WORKDIR}"

echo "========================================"
echo "WASP2-Rust Pipeline (chr22 subset)"
echo "Timestamp: $(date)"
echo "========================================"

# Python script for Rust module calls
python3 << 'PYTHON_SCRIPT'
import sys
import time
import json
sys.path.insert(0, '/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp')

import wasp2_rust
import subprocess

WORKDIR = "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison"
RUST_DIR = f"{WORKDIR}/rust"
INPUT_BAM = f"{WORKDIR}/input/chr22_input.bam"
INPUT_VCF = f"{WORKDIR}/input/chr22_variants.vcf.gz"
REFERENCE = f"{WORKDIR}/input/reference.fa"

results = {"pipeline": "wasp2_rust", "steps": {}}

print("\n" + "="*50)
print("STEP 1: VCF to BED")
print("="*50)
start = time.time()
bed_file = f"{RUST_DIR}/step1_vcf_to_bed/variants.bed"
wasp2_rust.vcf_to_bed(INPUT_VCF, bed_file, "NA12878")
elapsed = time.time() - start
results["steps"]["step1_vcf_to_bed"] = {"time_s": elapsed}

# Count variants
with open(bed_file) as f:
    variant_count = sum(1 for _ in f)
results["steps"]["step1_vcf_to_bed"]["variant_count"] = variant_count
print(f"  Variants in BED: {variant_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 2: Filter BAM by Variants")
print("="*50)
start = time.time()
to_remap_bam = f"{RUST_DIR}/step2_filter_bam/to_remap.bam"
keep_bam = f"{RUST_DIR}/step2_filter_bam/keep.bam"
wasp2_rust.filter_bam_by_variants(INPUT_BAM, bed_file, to_remap_bam, keep_bam)
elapsed = time.time() - start
results["steps"]["step2_filter_bam"] = {"time_s": elapsed}

# Count reads and extract read names
to_remap_count = int(subprocess.check_output(f"samtools view -c {to_remap_bam}", shell=True).decode().strip())
keep_count = int(subprocess.check_output(f"samtools view -c {keep_bam}", shell=True).decode().strip())
results["steps"]["step2_filter_bam"]["to_remap_reads"] = to_remap_count
results["steps"]["step2_filter_bam"]["keep_reads"] = keep_count

# Extract unique read names for comparison
subprocess.run(f"samtools view {to_remap_bam} | cut -f1 | sort -u > {RUST_DIR}/step2_filter_bam/read_names.txt", shell=True)
unique_names = int(subprocess.check_output(f"wc -l < {RUST_DIR}/step2_filter_bam/read_names.txt", shell=True).decode().strip())
results["steps"]["step2_filter_bam"]["unique_read_names"] = unique_names

print(f"  To remap reads: {to_remap_count}")
print(f"  Keep reads: {keep_count}")
print(f"  Unique read names to remap: {unique_names}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 3: Make Remap Reads (unified)")
print("="*50)
start = time.time()
remap_fastq = f"{RUST_DIR}/step3_make_reads/to_remap.fq.gz"
wasp2_rust.unified_make_reads(to_remap_bam, bed_file, remap_fastq, "NA12878")
elapsed = time.time() - start
results["steps"]["step3_make_reads"] = {"time_s": elapsed}

# Count FASTQ entries
fastq_count = int(subprocess.check_output(f"zcat {remap_fastq} | grep -c '^@'", shell=True).decode().strip())
results["steps"]["step3_make_reads"]["fastq_entries"] = fastq_count

# Extract WASP names for analysis
subprocess.run(f"zcat {remap_fastq} | grep '^@' | sed 's/^@//' | cut -d' ' -f1 > {RUST_DIR}/step3_make_reads/wasp_names.txt", shell=True)

print(f"  FASTQ entries: {fastq_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 4: BWA-MEM Remap")
print("="*50)
start = time.time()
remapped_bam = f"{RUST_DIR}/step4_remap/remapped.bam"
subprocess.run(f"bwa mem -t 4 -p {REFERENCE} {remap_fastq} 2>/dev/null | samtools sort -@ 4 -o {remapped_bam}", shell=True, check=True)
subprocess.run(f"samtools index {remapped_bam}", shell=True, check=True)
elapsed = time.time() - start
results["steps"]["step4_remap"] = {"time_s": elapsed}

remapped_count = int(subprocess.check_output(f"samtools view -c {remapped_bam}", shell=True).decode().strip())
results["steps"]["step4_remap"]["remapped_reads"] = remapped_count
print(f"  Remapped reads: {remapped_count}")
print(f"  Time: {elapsed:.2f}s")

print("\n" + "="*50)
print("STEP 5: Filter Remapped")
print("="*50)
start = time.time()
keep_remapped = f"{RUST_DIR}/step5_filter_remapped/keep.bam"
removed_bam = f"{RUST_DIR}/step5_filter_remapped/removed.bam"
wasp2_rust.filter_remapped_reads(remapped_bam, keep_remapped, removed_bam)
elapsed = time.time() - start
results["steps"]["step5_filter_remapped"] = {"time_s": elapsed}

kept_count = int(subprocess.check_output(f"samtools view -c {keep_remapped}", shell=True).decode().strip())
removed_count = int(subprocess.check_output(f"samtools view -c {removed_bam}", shell=True).decode().strip())
results["steps"]["step5_filter_remapped"]["kept_reads"] = kept_count
results["steps"]["step5_filter_remapped"]["removed_reads"] = removed_count

# Extract kept/removed read names
subprocess.run(f"samtools view {keep_remapped} | cut -f1 | sort -u > {RUST_DIR}/step5_filter_remapped/kept_names.txt", shell=True)
subprocess.run(f"samtools view {removed_bam} | cut -f1 | sort -u > {RUST_DIR}/step5_filter_remapped/removed_names.txt", shell=True)

print(f"  Kept reads: {kept_count}")
print(f"  Removed reads: {removed_count}")
print(f"  Time: {elapsed:.2f}s")

# Save results
with open(f"{RUST_DIR}/pipeline_results.json", 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*50)
print("RUST PIPELINE COMPLETE")
print("="*50)
print(f"Results saved to: {RUST_DIR}/pipeline_results.json")

PYTHON_SCRIPT

echo ""
echo "Completed: $(date)"
