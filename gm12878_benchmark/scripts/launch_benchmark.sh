#!/bin/bash
# GM12878 Benchmark Launcher with Indels
# This script launches the benchmark in the background using nohup

WORKDIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp"
SCRIPT="${WORKDIR}/gm12878_benchmark/scripts/run_benchmark.sh"
LOGFILE="${WORKDIR}/gm12878_benchmark/logs/nohup.log"

cd ${WORKDIR}

echo "Starting GM12878 benchmark with indels at $(date)"
echo "Working directory: ${WORKDIR}"
echo "Log file: ${LOGFILE}"

nohup bash ${SCRIPT} > ${LOGFILE} 2>&1 &

PID=$!
echo "Benchmark launched with PID: ${PID}"
echo "Monitor progress with: tail -f ${LOGFILE}"
echo "Check if running with: ps -p ${PID}"
