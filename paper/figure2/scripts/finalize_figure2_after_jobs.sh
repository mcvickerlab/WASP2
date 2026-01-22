#!/bin/bash
# Finalize Figure 2 after upstream cluster jobs finish.
#
# Usage:
#   qsub -hold_jid <figure2_bench_job>,<wasp1_counts_job> paper/figure2/scripts/finalize_figure2_after_jobs.sh
#
# This will:
#  - regenerate merged tables (count_comparison, before_after)
#  - regenerate Figure 2 with Panel C summary bars (and WASP1 included if present)
#
#$ -N fig2_finalize
#$ -cwd
#$ -j y
#$ -V
#$ -o paper/figure2/logs/fig2_finalize.$JOB_ID.log
#$ -pe iblm 2
#$ -l h_vmem=4G
#$ -l h_rt=1:00:00

set -eo pipefail

REPO_ROOT="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
cd "$REPO_ROOT"

# Make this job robust to inherited login-shell environment variables (-V).
# In particular, avoid surprising PYTHONPATH entries and avoid BLAS oversubscription.
unset PYTHONPATH || true
unset PYTHONHOME || true
export PYTHONPATH="${REPO_ROOT}/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

echo "[$(date)] Regenerating Figure 2 inputs/plots..."

python -c "import numpy as np; print('numpy', np.__version__, np.__file__)"

python paper/figure2/scripts/generate_count_comparison.py --dataset hg00731
python paper/figure2/scripts/generate_before_after_counts.py --dataset hg00731
python paper/figure2/scripts/generate_figure2.py --dataset hg00731 --panel-c-style summary_bars

# Optional: if WASP1 counts exist, emit the WASP1 vs WASP2 filtering supplement.
if [ -f paper/figure2/data/hg00731/wasp1_counts.filtered.tsv ]; then
  python paper/figure2/scripts/generate_wasp1_vs_wasp2_filtering.py --dataset hg00731
fi

echo "[$(date)] Done."
