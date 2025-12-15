#!/bin/bash
# Prepare HG00731 benchmark inputs (VCF + FASTQs) for STAR/WASP runs.
#
# - Downloads ERR1050079 FASTQs from ENA if missing
# - Builds `HG00731_het_only_chr.vcf.gz` (het-only, chr-prefixed) from 1000G
#
# This is intentionally separate from the benchmark scripts so failures are
# caught early and reruns don't repeatedly download large inputs.
#
#$ -N hg00731_data_prep
#$ -V
#$ -pe iblm 1
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

# NOTE: Avoid `set -u` here; conda activation scripts may reference unset vars.
set -eo pipefail

WASP2_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
DATA_DIR="${WASP2_DIR}/benchmarking/star_wasp_comparison/data"

mkdir -p "${DATA_DIR}"

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

echo "== HG00731 data prep =="
echo "Date: $(date)"
echo "Data dir: ${DATA_DIR}"

# -----------------------------------------------------------------------------
# FASTQs (ENA)
# -----------------------------------------------------------------------------
R1="${DATA_DIR}/ERR1050079_1.fastq.gz"
R2="${DATA_DIR}/ERR1050079_2.fastq.gz"
ENA_BASE="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/009/ERR1050079"

if [ ! -s "${R1}" ]; then
  echo "Downloading ${R1}..."
  wget -c --progress=dot:giga -O "${R1}" "${ENA_BASE}/ERR1050079_1.fastq.gz"
else
  echo "FASTQ present: ${R1}"
fi

if [ ! -s "${R2}" ]; then
  echo "Downloading ${R2}..."
  wget -c --progress=dot:giga -O "${R2}" "${ENA_BASE}/ERR1050079_2.fastq.gz"
else
  echo "FASTQ present: ${R2}"
fi

ls -lh "${R1}" "${R2}"

# -----------------------------------------------------------------------------
# VCF (1000G -> het-only -> chr-prefixed)
# -----------------------------------------------------------------------------
VCF_OUT="${DATA_DIR}/HG00731_het_only_chr.vcf.gz"
VCF_TBI="${VCF_OUT}.tbi"

if [ -s "${VCF_OUT}" ] && [ -s "${VCF_TBI}" ]; then
  echo "VCF present: ${VCF_OUT}"
else
  echo "Building ${VCF_OUT} (this can take a while)..."
  bash "${WASP2_DIR}/benchmarking/star_wasp_comparison/scripts/download_and_extract_vcf.sh"
fi

ls -lh "${VCF_OUT}" "${VCF_TBI}"
echo "Done."
