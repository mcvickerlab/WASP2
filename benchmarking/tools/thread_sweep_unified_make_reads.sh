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
Usage: $(basename "$0") --bam PATH --bed PATH --out-dir PATH [options]

Runs a thread/compression sweep for the Rust unified make-reads pipeline on an existing BAM+BED.
Each run spawns a fresh Python process to avoid cross-run state.

Required:
  --bam PATH              Coordinate-sorted, indexed BAM
  --bed PATH              Variant BED (vcf_to_bed output)
  --out-dir PATH          Output directory for TSV + per-run stats JSON

Optional:
  --label STR             Label for output row prefix (default: dataset basename)
  --threads LIST          Comma list (default: 1,2,4,8)
  --compression LIST      Comma list (default: 1,2,4)
  --compress-output 0/1   Whether to gzip FASTQs (default: 1)
  --indel-mode 0/1        Enable indel_mode (default: 0)
  --max-indel-size N      Only used when --indel-mode 1 (default: 50)
  --max-seqs N            (default: 64)
  --enable-timing 0/1     Export WASP2_TIMING=1 (adds per-read timing overhead; default: 0)

Outputs:
  - \${out-dir}/thread_sweep.tsv
  - \${out-dir}/runs/<params>/unified_stats.json
EOF
}

INPUT_BAM=""
INPUT_BED=""
OUT_DIR=""
LABEL=""
THREAD_LIST="1,2,4,8"
COMP_LIST="1,2,4"
COMPRESS_OUTPUT=1
INDEL_MODE=0
MAX_INDEL_SIZE=50
MAX_SEQS=64
ENABLE_TIMING=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) INPUT_BAM="$2"; shift 2 ;;
    --bed) INPUT_BED="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --threads) THREAD_LIST="$2"; shift 2 ;;
    --compression) COMP_LIST="$2"; shift 2 ;;
    --compress-output) COMPRESS_OUTPUT="$2"; shift 2 ;;
    --indel-mode) INDEL_MODE="$2"; shift 2 ;;
    --max-indel-size) MAX_INDEL_SIZE="$2"; shift 2 ;;
    --max-seqs) MAX_SEQS="$2"; shift 2 ;;
    --enable-timing) ENABLE_TIMING="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${INPUT_BAM}" || -z "${INPUT_BED}" || -z "${OUT_DIR}" ]]; then
  usage
  exit 2
fi

if [[ ! -f "${INPUT_BAM}" ]]; then
  echo "ERROR: --bam not found: ${INPUT_BAM}" >&2
  exit 2
fi
if [[ ! -f "${INPUT_BAM}.bai" && ! -f "${INPUT_BAM%.bam}.bai" ]]; then
  echo "ERROR: BAM index (.bai) not found for: ${INPUT_BAM}" >&2
  exit 2
fi
if [[ ! -f "${INPUT_BED}" ]]; then
  echo "ERROR: --bed not found: ${INPUT_BED}" >&2
  exit 2
fi

mkdir -p "${OUT_DIR}/runs"

if [[ -z "${LABEL}" ]]; then
  LABEL="$(basename "${INPUT_BAM}")"
fi

TSV="${OUT_DIR}/thread_sweep.tsv"
if [[ ! -f "${TSV}" ]]; then
  printf "label\tthreads\tcompression_threads\tcompress_output\tindel_mode\tmax_indel_size\tmax_seqs\twall_s\tpairs_processed\thaplotypes_written\ttree_build_ms\tbam_stream_ms\toverlap_query_ms\tpair_process_ms\tsend_ms\twriter_thread_ms\n" > "${TSV}"
fi

IFS=',' read -r -a THREADS_ARR <<< "${THREAD_LIST}"
IFS=',' read -r -a COMP_ARR <<< "${COMP_LIST}"

echo "Output: ${OUT_DIR}"
echo "BAM: ${INPUT_BAM}"
echo "BED: ${INPUT_BED}"
echo "Label: ${LABEL}"
echo "Threads: ${THREAD_LIST}"
echo "Compression threads: ${COMP_LIST}"
echo "compress_output: ${COMPRESS_OUTPUT}"
echo "indel_mode: ${INDEL_MODE}"
echo "max_indel_size: ${MAX_INDEL_SIZE}"
echo "max_seqs: ${MAX_SEQS}"
echo "enable_timing: ${ENABLE_TIMING}"
echo ""

for t in "${THREADS_ARR[@]}"; do
  for c in "${COMP_ARR[@]}"; do
    RUN_DIR="${OUT_DIR}/runs/t${t}_c${c}_gz${COMPRESS_OUTPUT}_indel${INDEL_MODE}"
    mkdir -p "${RUN_DIR}"

    OUT_R1="${RUN_DIR}/remap_r1.fq.gz"
    OUT_R2="${RUN_DIR}/remap_r2.fq.gz"
    if [[ "${COMPRESS_OUTPUT}" != "1" ]]; then
      OUT_R1="${RUN_DIR}/remap_r1.fq"
      OUT_R2="${RUN_DIR}/remap_r2.fq"
    fi

    echo "RUN threads=${t} compression_threads=${c} compress_output=${COMPRESS_OUTPUT} indel_mode=${INDEL_MODE}"

    # NOTE: WASP2_TIMING adds overhead (per-read Instant::now). Keep off for runtime sweeps.
    if [[ "${ENABLE_TIMING}" == "1" ]]; then
      export WASP2_TIMING=1
    else
      unset WASP2_TIMING || true
    fi

    python - <<PY
import json
import time
from pathlib import Path

from wasp2_rust import unified_make_reads_parallel_py

bam = "${INPUT_BAM}"
bed = "${INPUT_BED}"
out_r1 = "${OUT_R1}"
out_r2 = "${OUT_R2}"

t0 = time.time()
stats = unified_make_reads_parallel_py(
    bam,
    bed,
    out_r1,
    out_r2,
    max_seqs=int(${MAX_SEQS}),
    threads=int(${t}),
    compression_threads=int(${c}),
    compress_output=bool(int(${COMPRESS_OUTPUT})),
    indel_mode=bool(int(${INDEL_MODE})),
    max_indel_size=int(${MAX_INDEL_SIZE}),
)
wall_s = time.time() - t0
stats = dict(stats)
stats["wall_s"] = wall_s

out_json = Path("${RUN_DIR}") / "unified_stats.json"
out_json.write_text(json.dumps(stats, indent=2))
print(json.dumps({"pairs_processed": stats.get("pairs_processed"), "haplotypes_written": stats.get("haplotypes_written"), "wall_s": wall_s}))
PY

    # Append TSV row from per-run JSON
    python - <<PY >> "${TSV}"
import json

d = json.load(open("${RUN_DIR}/unified_stats.json"))
row = [
  "${LABEL}",
  "${t}",
  "${c}",
  str(int(${COMPRESS_OUTPUT})),
  str(int(${INDEL_MODE})),
  str(int(${MAX_INDEL_SIZE})),
  str(int(${MAX_SEQS})),
  str(d.get("wall_s", "NA")),
  str(d.get("pairs_processed", "NA")),
  str(d.get("haplotypes_written", "NA")),
  str(d.get("tree_build_ms", "NA")),
  str(d.get("bam_stream_ms", "NA")),
  str(d.get("overlap_query_ms", "NA")),
  str(d.get("pair_process_ms", "NA")),
  str(d.get("send_ms", "NA")),
  str(d.get("writer_thread_ms", "NA")),
]
print("\t".join(row))
PY

    # Remove potentially large FASTQs/sidecar; keep only stats JSON.
    rm -f "${OUT_R1}" "${OUT_R2}" "${OUT_R1}.expected_positions.tsv" 2>/dev/null || true
  done
done

echo ""
echo "Done. Results: ${TSV}"
