#!/bin/bash
# Wait for an SGE array job to finish, then rebuild the INDEL scaling COMBINED.tsv.
#
# Usage:
#   benchmarking/wait_and_combine_indel_scaling.sh [JOB_ID] [SLEEP_SECONDS]
#
# Defaults:
#   JOB_ID=8724785 (current indel scaling array)
#   SLEEP_SECONDS=60

set -euo pipefail

JOB_ID="${1:-8724785}"
SLEEP_SECONDS="${2:-60}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "Waiting for job ${JOB_ID} to finish..."
while qstat 2>/dev/null | awk 'NR>2{print $1}' | grep -q "^${JOB_ID}$"; do
  echo "$(date) job ${JOB_ID} still running/queued; sleeping ${SLEEP_SECONDS}s"
  sleep "${SLEEP_SECONDS}"
done

echo "$(date) job ${JOB_ID} finished; combining results"
python "${SCRIPT_DIR}/combine_unified_indel_scaling.py" \
  --out "${ROOT_DIR}/benchmarking/results/wasp2_unified_indel_scaling_COMBINED.tsv"

echo "Done."

