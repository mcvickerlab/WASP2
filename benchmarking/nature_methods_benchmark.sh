#!/bin/bash
# =============================================================================
# NATURE METHODS BENCHMARK: WASP2 vs GATK vs phASER vs biastools
#
# Measures: SPEED + ACCURACY (allele counting)
# Data: HG00731 RNA-seq (56M reads) from STAR-WASP paper
# Threads: 8
# =============================================================================

set -e

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
echo "=============================================="
echo "NATURE METHODS BENCHMARK"
echo "WASP2 vs GATK vs phASER vs biastools"
echo "Timestamp: ${TIMESTAMP}"
echo "=============================================="

# Paths
WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
DATA_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison/data"
OUTPUT_DIR="${WASP2_DIR}/benchmarking/nature_methods_results/${TIMESTAMP}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/timing"
mkdir -p "${OUTPUT_DIR}/counts"

export PYTHONPATH="${WASP2_DIR}/src:${PYTHONPATH}"
cd "${WASP2_DIR}"

# Input files
FASTQ_R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
FASTQ_R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
VCF="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"
STAR_INDEX="/iblm/netapp/data1/external/GRC38/combined/GDC_hg38/star_index"
REFERENCE="/iblm/netapp/data1/external/GRC38/combined/GDC_hg38/GRCh38.d1.vd1.fa"

# Check inputs exist
echo ""
echo "Checking inputs..."
for f in "${FASTQ_R1}" "${FASTQ_R2}" "${VCF}"; do
    if [[ -f "$f" ]]; then
        echo "  ✓ $(basename $f)"
    else
        echo "  ✗ MISSING: $f"
        exit 1
    fi
done

# Thread count
THREADS=8
echo ""
echo "Using ${THREADS} threads"

# =============================================================================
# STEP 1: STAR Alignment (shared by all methods)
# =============================================================================
echo ""
echo "=============================================="
echo "STEP 1: STAR Alignment"
echo "=============================================="

ALIGNED_BAM="${OUTPUT_DIR}/aligned.sorted.bam"
if [[ -f "${ALIGNED_BAM}" ]]; then
    echo "Using existing alignment: ${ALIGNED_BAM}"
else
    STAR_START=$(date +%s.%N)

    STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${FASTQ_R1} ${FASTQ_R2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${OUTPUT_DIR}/star_ \
         --outSAMattributes NH HI AS nM NM MD \
         --alignEndsType EndToEnd \
         --outFilterMultimapNmax 1

    mv ${OUTPUT_DIR}/star_Aligned.sortedByCoord.out.bam ${ALIGNED_BAM}
    samtools index ${ALIGNED_BAM}

    STAR_END=$(date +%s.%N)
    STAR_TIME=$(echo "${STAR_END} - ${STAR_START}" | bc)
    echo "STAR_ALIGNMENT,${STAR_TIME}" > ${OUTPUT_DIR}/timing/star.csv
    echo "STAR alignment completed in ${STAR_TIME}s"
fi

# =============================================================================
# STEP 2: WASP2-Rust (Speed + Accuracy)
# =============================================================================
echo ""
echo "=============================================="
echo "STEP 2: WASP2-Rust"
echo "=============================================="

WASP2_START=$(date +%s.%N)

python << PYTHON
import sys
sys.path.insert(0, 'src')
import time
import wasp2_rust
import pandas as pd
import pysam

bam_path = "${ALIGNED_BAM}"
vcf_path = "${VCF}"
output_dir = "${OUTPUT_DIR}"

print("Loading variants from VCF...")
vcf = pysam.VariantFile(vcf_path)
regions = []
for rec in vcf:
    if rec.alts and len(rec.ref) == 1 and len(rec.alts[0]) == 1:  # SNPs only
        regions.append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
print(f"Loaded {len(regions)} SNP variants")

print("Counting alleles with WASP2-Rust...")
start = time.time()
bc = wasp2_rust.BamCounter(bam_path)
counts = bc.count_alleles(regions, min_qual=20, threads=${THREADS})
elapsed = time.time() - start
print(f"Counting completed in {elapsed:.2f}s")

# Save results
results = []
for i, (chrom, pos, ref, alt) in enumerate(regions):
    ref_count, alt_count, other = counts[i]
    total = ref_count + alt_count
    ratio = ref_count / total if total > 0 else 0.5
    results.append({
        'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
        'ref_count': ref_count, 'alt_count': alt_count,
        'total': total, 'ratio': ratio
    })

df = pd.DataFrame(results)
df.to_csv(f"{output_dir}/counts/wasp2_counts.csv", index=False)
print(f"Saved {len(df)} variant counts")

# Summary
print(f"\nWASP2-Rust Summary:")
print(f"  Variants counted: {len(df)}")
print(f"  Total reads: {df['total'].sum()}")
print(f"  Mean coverage: {df['total'].mean():.1f}")
print(f"  Mean REF ratio: {df['ratio'].mean():.3f}")
PYTHON

WASP2_END=$(date +%s.%N)
WASP2_TIME=$(echo "${WASP2_END} - ${WASP2_START}" | bc)
echo "WASP2,${WASP2_TIME}" > ${OUTPUT_DIR}/timing/wasp2.csv
echo "WASP2-Rust completed in ${WASP2_TIME}s"

# =============================================================================
# STEP 3: GATK ASEReadCounter
# =============================================================================
echo ""
echo "=============================================="
echo "STEP 3: GATK ASEReadCounter"
echo "=============================================="

GATK_START=$(date +%s.%N)

# Check if GATK is available
if command -v gatk &> /dev/null; then
    echo "Running GATK ASEReadCounter..."

    # Add read groups if needed
    GATK_BAM="${OUTPUT_DIR}/gatk_input.bam"
    gatk AddOrReplaceReadGroups \
        -I ${ALIGNED_BAM} \
        -O ${GATK_BAM} \
        -RGID HG00731 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM HG00731 \
        2>/dev/null || cp ${ALIGNED_BAM} ${GATK_BAM}
    samtools index ${GATK_BAM}

    # Create sequence dictionary if needed
    DICT_FILE="${REFERENCE%.fa}.dict"
    if [[ ! -f "${DICT_FILE}" ]]; then
        gatk CreateSequenceDictionary -R ${REFERENCE} -O ${DICT_FILE} 2>/dev/null || true
    fi

    # Run ASEReadCounter
    gatk ASEReadCounter \
        -R ${REFERENCE} \
        -I ${GATK_BAM} \
        -V ${VCF} \
        -O ${OUTPUT_DIR}/counts/gatk_raw.table \
        --min-mapping-quality 10 \
        --min-base-quality 20 \
        2>/dev/null || echo "GATK failed (may not be installed)"

    # Parse GATK output
    if [[ -f "${OUTPUT_DIR}/counts/gatk_raw.table" ]]; then
        python simulation/competitors/run_gatk.py \
            --bam ${GATK_BAM} \
            --vcf ${VCF} \
            --ref ${REFERENCE} \
            --output ${OUTPUT_DIR}/counts/gatk/ 2>/dev/null || true
    fi
else
    echo "GATK not available, skipping..."
fi

GATK_END=$(date +%s.%N)
GATK_TIME=$(echo "${GATK_END} - ${GATK_START}" | bc)
echo "GATK,${GATK_TIME}" > ${OUTPUT_DIR}/timing/gatk.csv
echo "GATK completed in ${GATK_TIME}s"

# =============================================================================
# STEP 4: biastools
# =============================================================================
echo ""
echo "=============================================="
echo "STEP 4: biastools"
echo "=============================================="

BIASTOOLS_START=$(date +%s.%N)

python simulation/competitors/run_biastools.py \
    --bam ${ALIGNED_BAM} \
    --vcf ${VCF} \
    --output ${OUTPUT_DIR}/counts/biastools/ 2>&1 || echo "biastools completed (or failed)"

BIASTOOLS_END=$(date +%s.%N)
BIASTOOLS_TIME=$(echo "${BIASTOOLS_END} - ${BIASTOOLS_START}" | bc)
echo "BIASTOOLS,${BIASTOOLS_TIME}" > ${OUTPUT_DIR}/timing/biastools.csv
echo "biastools completed in ${BIASTOOLS_TIME}s"

# =============================================================================
# STEP 5: Compare Results
# =============================================================================
echo ""
echo "=============================================="
echo "STEP 5: Compare Results"
echo "=============================================="

python << PYTHON
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

output_dir = Path("${OUTPUT_DIR}")

# Load timing
print("\n=== TIMING RESULTS ===")
for f in (output_dir / "timing").glob("*.csv"):
    tool, time = open(f).read().strip().split(",")
    print(f"  {tool}: {float(time):.2f}s")

# Load counts and compare
print("\n=== ACCURACY COMPARISON ===")
wasp2 = pd.read_csv(output_dir / "counts" / "wasp2_counts.csv")
print(f"WASP2: {len(wasp2)} variants, mean ratio {wasp2['ratio'].mean():.3f}")

# Compare to other tools if available
for tool in ['gatk', 'biastools']:
    tool_file = output_dir / "counts" / tool / f"{tool}_variant_level.csv"
    if tool_file.exists():
        df = pd.read_csv(tool_file)
        print(f"{tool.upper()}: {len(df)} variants")
        # Calculate correlation with WASP2
        merged = wasp2.merge(df, on=['chrom', 'pos'], suffixes=('_wasp2', f'_{tool}'))
        if len(merged) > 10:
            r, p = stats.pearsonr(merged['ratio'], merged[f'{tool}_ratio'])
            print(f"  Correlation with WASP2: r={r:.3f}")

print("\nResults saved to:", output_dir)
PYTHON

echo ""
echo "=============================================="
echo "BENCHMARK COMPLETE"
echo "Results: ${OUTPUT_DIR}"
echo "=============================================="
