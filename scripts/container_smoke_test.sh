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
echo "[1/5] Version checks..."

for cmd in wasp2-count wasp2-map wasp2-analyze; do
    if $cmd --version > /dev/null 2>&1; then
        echo "  PASS: $cmd --version"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $cmd --version"
        FAIL=$((FAIL + 1))
    fi
done
python - <<'PY'
import os
from importlib.metadata import version
import wasp2
import wasp2_rust

installed = version("wasp2")
assert wasp2.__version__ == wasp2_rust.__version__ == installed
expected = os.environ.get("EXPECTED_WASP2_VERSION", installed)
assert installed == expected, (installed, expected)
print(f"  PASS: Python, Rust and distribution versions agree: {installed}")
PY
PASS=$((PASS + 1))
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Python imports
# ─────────────────────────────────────────────────────────────────────────────
echo "[2/5] Python imports..."

if python -c "import wasp2_rust; print('Rust extension OK')" 2>/dev/null; then
    echo "  PASS: wasp2_rust import"
    PASS=$((PASS + 1))
else
    echo "  FAIL: wasp2_rust import"
    FAIL=$((FAIL + 1))
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
# Test 3: Count variants with the bundled synthetic BAM/VCF
# ─────────────────────────────────────────────────────────────────────────────
echo "[3/5] Count variants with bundled BAM/VCF..."

if [[ -f "$DATA_DIR/sample1.bam" && -f "$DATA_DIR/variants.vcf.gz" ]]; then
    if wasp2-count count-variants \
        "$DATA_DIR/sample1.bam" \
        "$DATA_DIR/variants.vcf.gz" \
        --samples sample1 \
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
echo "[4/5] External tool availability..."

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

# Exercise the original ATAC failure: reordered columns, nonnumeric region IDs,
# and both phase orientations. Compare phased and unphased CLI runs to catch
# silently ignored --phased as well as positional parsing and ignored grouping.
echo "[5/5] ATAC region analysis..."
python - "$TMP_DIR/regions.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], "w") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["ref", "alt", "region", "chrom", "start", "pos", "GT", "ref_count", "alt_count"])
    for pos, gt, ref, alt in [(101, "0|1", 18, 2), (201, "1|0", 3, 17),
                              (301, "0|1", 22, 18), (401, "1|0", 25, 15)]:
        region = "chr1:100-250" if pos < 300 else "chr1:300-450"
        writer.writerow(["A", "G", region, "chr1", pos - 1, pos, gt, ref, alt])
PY
wasp2-analyze find-imbalance "$TMP_DIR/regions.tsv" --phased \
    --region_col region --out_file "$TMP_DIR/phased.tsv"
wasp2-analyze find-imbalance "$TMP_DIR/regions.tsv" \
    --region_col region --out_file "$TMP_DIR/unphased.tsv"
python - "$TMP_DIR/phased.tsv" "$TMP_DIR/unphased.tsv" <<'PY'
import csv
import math
import sys

def read(path):
    with open(path) as handle:
        return {row["region"]: row for row in csv.DictReader(handle, delimiter="\t")}

phased, unphased = map(read, sys.argv[1:])
assert set(phased) == set(unphased) == {"chr1:100-250", "chr1:300-450"}
for region, counts in [("chr1:100-250", (21, 19)), ("chr1:300-450", (47, 33))]:
    for results in (phased, unphased):
        row = results[region]
        assert (int(row["ref_count"]), int(row["alt_count"])) == counts
        assert int(row["snp_count"]) == 2
        assert all(math.isfinite(float(row[key])) for key in ["mu", "null_ll", "alt_ll", "pval", "fdr_pval"])
assert any(abs(float(phased[k]["alt_ll"]) - float(unphased[k]["alt_ll"])) > 1e-6 for k in phased)
print("  PASS: region grouping, raw counts and phased likelihood through CLI")
PY
PASS=$((PASS + 1))
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
