#!/bin/bash
# =============================================================================
# WASP2 Container Smoke Test
# =============================================================================
# Validates WASP2 tools inside a Docker/Singularity container.
# Designed to be bundled at /opt/wasp2/scripts/container_smoke_test.sh
# with test data at /opt/wasp2/test-data/ (copied from tests/shared_data/).
#
# When run outside a container, uses tests/shared_data/ directly.
#
# Usage (container):
#   docker run wasp2:test /opt/wasp2/scripts/container_smoke_test.sh
#   singularity exec wasp2.sif /opt/wasp2/scripts/container_smoke_test.sh
#
# Usage (local):
#   bash scripts/container_smoke_test.sh
#
# Exit codes:
#   0 = all tests passed
#   1 = one or more tests failed
# =============================================================================

set -euo pipefail

# Detect data directory: container path or repo path
if [[ -d "/opt/wasp2/test-data" ]]; then
    DATA_DIR="/opt/wasp2/test-data"
elif [[ -d "tests/shared_data" ]]; then
    DATA_DIR="tests/shared_data"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DATA_DIR="$(cd "$SCRIPT_DIR/../tests/shared_data" 2>/dev/null && pwd)" || {
        echo "ERROR: Cannot find test data directory"
        exit 1
    }
fi

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/wasp2_container_smoke_XXXXXX")
PASS=0
FAIL=0

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

echo "==================================================================="
echo " WASP2 Container Smoke Test"
echo "==================================================================="
echo "Data dir: $DATA_DIR"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: CLI binaries exist and print version
# ─────────────────────────────────────────────────────────────────────────────
echo "[1/4] Version checks..."

for cmd in wasp2-count wasp2-map wasp2-analyze; do
    if $cmd --version > /dev/null 2>&1; then
        echo "  PASS: $cmd --version"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $cmd --version"
        FAIL=$((FAIL + 1))
    fi
done
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Python imports
# ─────────────────────────────────────────────────────────────────────────────
echo "[2/4] Python imports..."

if python -c "import wasp2_rust; print('Rust extension OK')" 2>/dev/null; then
    echo "  PASS: wasp2_rust import"
    PASS=$((PASS + 1))
else
    echo "  INFO: wasp2_rust not available (Python-only mode)"
fi

if python -c "from counting.run_counting import run_count_variants; print('counting OK')" 2>/dev/null; then
    echo "  PASS: counting module import"
    PASS=$((PASS + 1))
else
    echo "  FAIL: counting module import"
    FAIL=$((FAIL + 1))
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Count variants with real data
# ─────────────────────────────────────────────────────────────────────────────
echo "[3/4] Count variants with real data..."

if [[ -f "$DATA_DIR/sample1.bam" && -f "$DATA_DIR/variants.vcf.gz" ]]; then
    if wasp2-count count-variants \
        "$DATA_DIR/sample1.bam" \
        "$DATA_DIR/variants.vcf.gz" \
        --samples SAMPLE1 \
        --out "$TMP_DIR/counts.tsv" \
        2>/dev/null; then

        if [[ -s "$TMP_DIR/counts.tsv" ]]; then
            ROWS=$(wc -l < "$TMP_DIR/counts.tsv")
            echo "  PASS: count-variants produced $ROWS lines"
            PASS=$((PASS + 1))
        else
            echo "  FAIL: count-variants output is empty"
            FAIL=$((FAIL + 1))
        fi
    else
        echo "  FAIL: count-variants exited with error"
        FAIL=$((FAIL + 1))
    fi
else
    echo "  SKIP: Test data not available"
fi
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: External tools (samtools, bcftools, bedtools)
# ─────────────────────────────────────────────────────────────────────────────
echo "[4/4] External tool availability..."

for tool in samtools bcftools bedtools; do
    if command -v $tool > /dev/null 2>&1; then
        echo "  PASS: $tool available"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $tool not found"
        FAIL=$((FAIL + 1))
    fi
done
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
TOTAL=$((PASS + FAIL))
echo "==================================================================="
echo " Results: $PASS/$TOTAL passed"
echo "==================================================================="

if [[ $FAIL -gt 0 ]]; then
    echo " WARNING: $FAIL test(s) failed"
    exit 1
else
    echo " All container smoke tests passed!"
    exit 0
fi
