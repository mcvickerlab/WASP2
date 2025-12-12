#!/bin/bash
# SGE wrapper to submit WASP1 benchmark
# Usage: qsub sge_wasp1_benchmark.sh

cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking

# Ensure logs directory exists
mkdir -p logs

# Submit the WASP1 benchmark
qsub star_wasp_comparison/scripts/run_wasp1_benchmark.sh

echo "WASP1 benchmark submitted. Check logs/ for output."
