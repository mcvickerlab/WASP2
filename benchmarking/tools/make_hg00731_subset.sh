#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WASP2_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Conda activation scripts may reference unset vars (e.g. target_platform).
set +u
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2
set -u

usage() {
  cat <<EOF
Usage: $(basename "$0") [--bam PATH] [--bed PATH] [--out-dir PATH] [--threads N] [--seed N]

Creates reproducible HG00731 RNA-seq subsets for profiling/thread sweeps:
  - chr15-22 subset BAM + index
  - 2% and 10% subsamples (pair-consistent; samtools --subsample)
  - chr15-22 subset of the provided BED

Defaults:
  --bam     ${WASP2_DIR}/benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-10_19-23-02/original.bam
  --bed     ${WASP2_DIR}/benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-10_19-23-02/variants.bed
  --out-dir \${TMPDIR:-/tmp}/wasp2_hg00731_subset_<timestamp>
  --threads 8
  --seed    42
EOF
}

INPUT_BAM="${WASP2_DIR}/benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-10_19-23-02/original.bam"
INPUT_BED="${WASP2_DIR}/benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-10_19-23-02/variants.bed"
OUT_DIR=""
THREADS=8
SEED=42

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) INPUT_BAM="$2"; shift 2 ;;
    --bed) INPUT_BED="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${OUT_DIR}" ]]; then
  OUT_DIR="${TMPDIR:-/tmp}/wasp2_hg00731_subset_$(date +%Y%m%d_%H%M%S)"
fi

if [[ ! -f "${INPUT_BAM}" ]]; then
  echo "ERROR: --bam not found: ${INPUT_BAM}" >&2
  exit 2
fi
if [[ ! -f "${INPUT_BED}" ]]; then
  echo "ERROR: --bed not found: ${INPUT_BED}" >&2
  exit 2
fi

mkdir -p "${OUT_DIR}"

CHROMS=(chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

echo "Output: ${OUT_DIR}"
echo "Input BAM: ${INPUT_BAM}"
echo "Input BED: ${INPUT_BED}"
echo "Chroms: ${CHROMS[*]}"
echo "Threads: ${THREADS}"
echo "Seed: ${SEED}"

BASE_BAM="${OUT_DIR}/HG00731_chr15-22.bam"
echo ""
echo "[1/3] Extract chr15-22 BAM..."
samtools view -@ "${THREADS}" -b -o "${BASE_BAM}" "${INPUT_BAM}" "${CHROMS[@]}"
samtools index -@ "${THREADS}" "${BASE_BAM}"

echo ""
echo "[2/3] Subsample 2% and 10% (pair-consistent)..."
samtools view -@ "${THREADS}" -b --subsample 0.02 --subsample-seed "${SEED}" -o "${OUT_DIR}/HG00731_chr15-22_sub2pct.bam" "${BASE_BAM}"
samtools index -@ "${THREADS}" "${OUT_DIR}/HG00731_chr15-22_sub2pct.bam"

samtools view -@ "${THREADS}" -b --subsample 0.10 --subsample-seed "${SEED}" -o "${OUT_DIR}/HG00731_chr15-22_sub10pct.bam" "${BASE_BAM}"
samtools index -@ "${THREADS}" "${OUT_DIR}/HG00731_chr15-22_sub10pct.bam"

echo ""
echo "[3/3] Filter BED to chr15-22..."
OUT_BED="${OUT_DIR}/HG00731_chr15-22_het_only.bed"
awk -v chroms="${CHROMS[*]}" '
  BEGIN{
    n=split(chroms, a, " ");
    for(i=1;i<=n;i++) keep[a[i]]=1;
  }
  ($1 in keep){ print }
' "${INPUT_BED}" > "${OUT_BED}"

echo ""
echo "Done."
echo "BAMs:"
ls -lh "${OUT_DIR}"/HG00731_chr15-22*.bam
echo ""
echo "BED:"
ls -lh "${OUT_BED}"
