#!/bin/bash
#$ -N fig2_star_only
#$ -cwd
#$ -j y
#$ -V
#$ -o paper/figure2/logs/fig2_star_only.log
#$ -pe iblm 8
#$ -l h_vmem=32G
#$ -l h_rt=7:00:00
#
# Build an "original" (pre-WASP) HG00731 RNA-seq BAM for Figure 2 Panel C.
#
# IMPORTANT: Use the same STAR parameters as the FAIR WASP2-Rust benchmark
# (`benchmarking/star_wasp_comparison/scripts/run_wasp2_unified_FAIR.sh`)
# so "original vs WASP-filtered" comparisons reflect WASP filtering/remapping
# rather than differences in alignment settings.
#
# Output: `paper/figure2/data/hg00731/original.bam`

# NOTE: Avoid `set -u` here because conda's shell hooks are not nounset-safe.
set -eo pipefail

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT"

mkdir -p paper/figure2/data/hg00731 paper/figure2/logs

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

THREADS="${NSLOTS:-8}"
ulimit -n 10000

FASTQ_R1="$REPO_ROOT/benchmarking/star_wasp_comparison/data/ERR1050079_1.fastq.gz"
FASTQ_R2="$REPO_ROOT/benchmarking/star_wasp_comparison/data/ERR1050079_2.fastq.gz"
# Match the STAR index used by the existing RNA-seq benchmark pipelines.
STAR_INDEX="/iblm/netapp/data1/external/GRC38/combined/google_cloud/star_index"

OUT_PREFIX="$REPO_ROOT/paper/figure2/data/hg00731/star_only_"
OUT_BAM="$REPO_ROOT/paper/figure2/data/hg00731/original.bam"
STAR_BIN="STAR"

if [ ! -f "$FASTQ_R1" ] || [ ! -f "$FASTQ_R2" ]; then
  echo "ERROR: Missing FASTQs: $FASTQ_R1 or $FASTQ_R2" >&2
  exit 1
fi
if [ ! -d "$STAR_INDEX" ]; then
  echo "ERROR: Missing STAR index dir: $STAR_INDEX" >&2
  exit 1
fi
if ! command -v "$STAR_BIN" >/dev/null 2>&1; then
  echo "ERROR: STAR not found in PATH (WASP2_dev2)" >&2
  exit 1
fi
echo "Using STAR binary from PATH: $(command -v "$STAR_BIN")"
if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in PATH (WASP2_dev2)" >&2
  exit 1
fi

echo "Running STAR alignment (HG00731) with $THREADS threads (FAIR params)..."
echo "FASTQs: $FASTQ_R1 $FASTQ_R2"
echo "Index:  $STAR_INDEX"

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --limitBAMsortRAM 15000000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM NM MD jM jI \
  --alignEndsType EndToEnd \
  --outSAMunmapped Within \
  --outFilterMultimapNmax 1 \
  --readFilesIn "$FASTQ_R1" "$FASTQ_R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT_PREFIX"

echo "Finalizing BAM..."
mv "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" "$OUT_BAM"
samtools index -@ "$THREADS" "$OUT_BAM"

echo "Done:"
ls -lh "$OUT_BAM" "${OUT_BAM}.bai"
