#!/bin/bash
# Finalize Figure 1 inputs after reruns:
# - Combine unified scaling TSVs (SNV + INDEL) and copy into paper benchmarks
# - Copy latest RNA-seq WASP2-Rust JSONs into paper benchmarks
# - Regenerate `paper/figure1/plots/figure1.*`
#
# Intended usage:
#   qsub -hold_jid <snv_array>,<indel_array>,<rnaseq_snv_job>,<rnaseq_indel_job> benchmarking/qsub_finalize_figure1_reruns.sh
#
#$ -N fig1_finalize
#$ -V
#$ -pe iblm 1
#$ -l h_vmem=16G
#$ -j y
#$ -o /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/logs/
#$ -cwd

set -eo pipefail

ROOT_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
FIG_BENCH_DIR="${ROOT_DIR}/paper/figure1/data/benchmarks"

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

echo "== Figure 1 finalize =="
echo "Date: $(date)"
echo "ROOT_DIR: ${ROOT_DIR}"

cd "${ROOT_DIR}"

snv_scaling_dir="$(ls -d benchmarking/results/unified_scaling_* 2>/dev/null | sort | tail -1)"
indel_scaling_dir="$(ls -d benchmarking/results/unified_indel_scaling_* 2>/dev/null | sort | tail -1)"

if [ -z "${snv_scaling_dir}" ] || [ ! -f "${snv_scaling_dir}/wasp2_unified_scaling.tsv" ]; then
  echo "ERROR: missing unified scaling TSV under ${snv_scaling_dir}" >&2
  exit 1
fi
if [ -z "${indel_scaling_dir}" ] || [ ! -f "${indel_scaling_dir}/wasp2_unified_indel_scaling.tsv" ]; then
  echo "ERROR: missing unified INDEL scaling TSV under ${indel_scaling_dir}" >&2
  exit 1
fi

echo "SNV scaling dir:   ${snv_scaling_dir}"
echo "INDEL scaling dir: ${indel_scaling_dir}"

python benchmarking/combine_unified_scaling.py \
  --pattern "${snv_scaling_dir}/wasp2_unified_scaling.tsv" \
  --out benchmarking/results/wasp2_unified_scaling_COMBINED.tsv

python benchmarking/combine_unified_indel_scaling.py \
  --pattern "${indel_scaling_dir}/wasp2_unified_indel_scaling.tsv" \
  --out benchmarking/results/wasp2_unified_indel_scaling_COMBINED.tsv

mkdir -p "${FIG_BENCH_DIR}"
cp benchmarking/results/wasp2_unified_scaling_COMBINED.tsv "${FIG_BENCH_DIR}/panel_b_wasp2rust_snv_scaling.tsv"
cp benchmarking/results/wasp2_unified_indel_scaling_COMBINED.tsv "${FIG_BENCH_DIR}/panel_b_wasp2rust_indel_scaling.tsv"

echo "Updated Panel B benchmarks:"
ls -lh "${FIG_BENCH_DIR}/panel_b_wasp2rust_snv_scaling.tsv" "${FIG_BENCH_DIR}/panel_b_wasp2rust_indel_scaling.tsv"

rnaseq_dir_base="benchmarking/star_wasp_comparison/results"
rnaseq_snv_dir="$(ls -d ${rnaseq_dir_base}/wasp2rust_fair_* 2>/dev/null | sort | tail -1)"
rnaseq_indel_dir="$(ls -d ${rnaseq_dir_base}/wasp2rust_indel_rnaseq_v2_* 2>/dev/null | sort | tail -1)"

if [ -z "${rnaseq_snv_dir}" ] || [ ! -f "${rnaseq_snv_dir}/benchmark_results.json" ]; then
  echo "ERROR: missing RNA-seq SNV benchmark_results.json under ${rnaseq_snv_dir}" >&2
  exit 1
fi
if [ -z "${rnaseq_indel_dir}" ] || [ ! -f "${rnaseq_indel_dir}/benchmark_results.json" ]; then
  echo "ERROR: missing RNA-seq INDEL benchmark_results.json under ${rnaseq_indel_dir}" >&2
  exit 1
fi

cp "${rnaseq_snv_dir}/benchmark_results.json" "${FIG_BENCH_DIR}/panel_de_rnaseq_wasp2rust_snv.json"
cp "${rnaseq_indel_dir}/benchmark_results.json" "${FIG_BENCH_DIR}/panel_de_rnaseq_wasp2rust_indel.json"

echo "Updated Panel D/E benchmarks:"
ls -lh "${FIG_BENCH_DIR}/panel_de_rnaseq_wasp2rust_snv.json" "${FIG_BENCH_DIR}/panel_de_rnaseq_wasp2rust_indel.json"

echo "Regenerating Figure 1..."
python paper/figure1/scripts/generate_figure1.py

echo "Done. Outputs:"
ls -lh paper/figure1/plots/figure1.png paper/figure1/plots/figure1.pdf 2>/dev/null || true

