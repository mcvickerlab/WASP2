#!/bin/bash
# =============================================================================
# WASP2 Release Validation Script
# Run this BEFORE pushing to mcvickerlab/WASP2 upstream
# =============================================================================

set -e  # Exit on first error

echo "=============================================="
echo "WASP2 Release Validation"
echo "=============================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PASSED=0
FAILED=0

check_pass() {
    echo -e "${GREEN}✓ PASS${NC}: $1"
    ((PASSED++))
}

check_fail() {
    echo -e "${RED}✗ FAIL${NC}: $1"
    ((FAILED++))
}

check_warn() {
    echo -e "${YELLOW}⚠ WARN${NC}: $1"
}

# -----------------------------------------------------------------------------
# 1. Repository Size Check
# -----------------------------------------------------------------------------
echo "1. Checking repository size..."
REPO_SIZE=$(du -sm . --exclude=.git --exclude=rust/target 2>/dev/null | cut -f1)
if [ "$REPO_SIZE" -lt 200 ]; then
    check_pass "Repository size is ${REPO_SIZE}MB (expected <200MB for software-only)"
else
    check_fail "Repository size is ${REPO_SIZE}MB (too large - paper/benchmark files present?)"
fi

# -----------------------------------------------------------------------------
# 2. No Large Files Check
# -----------------------------------------------------------------------------
echo ""
echo "2. Checking for large files (>50MB)..."
LARGE_FILES=$(find . -type f -size +50M ! -path "./.git/*" ! -path "./rust/target/*" 2>/dev/null | wc -l)
if [ "$LARGE_FILES" -eq 0 ]; then
    check_pass "No large files found"
else
    check_fail "Found $LARGE_FILES files >50MB"
    find . -type f -size +50M ! -path "./.git/*" ! -path "./rust/target/*" 2>/dev/null
fi

# -----------------------------------------------------------------------------
# 3. Paper/Benchmark Directories Check (should NOT exist in release)
# -----------------------------------------------------------------------------
echo ""
echo "3. Checking that paper/benchmark directories are removed..."
UNWANTED_DIRS=("paper" "benchmarking" "results" "validation" "comparison_results")
FOUND_UNWANTED=0
for dir in "${UNWANTED_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        check_fail "Directory '$dir' exists (should be removed for release)"
        ((FOUND_UNWANTED++))
    fi
done
if [ "$FOUND_UNWANTED" -eq 0 ]; then
    check_pass "No paper/benchmark directories found"
fi

# -----------------------------------------------------------------------------
# 4. Required Directories Check
# -----------------------------------------------------------------------------
echo ""
echo "4. Checking required directories exist..."
REQUIRED_DIRS=("src" "rust" "tests" "docs")
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        check_pass "Directory '$dir' exists"
    else
        check_fail "Required directory '$dir' missing"
    fi
done

# -----------------------------------------------------------------------------
# 5. Required Files Check
# -----------------------------------------------------------------------------
echo ""
echo "5. Checking required files exist..."
REQUIRED_FILES=("pyproject.toml" "README.md" "LICENSE" "environment.yml")
for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        check_pass "File '$file' exists"
    else
        check_fail "Required file '$file' missing"
    fi
done

# -----------------------------------------------------------------------------
# 6. Python Tests
# -----------------------------------------------------------------------------
echo ""
echo "6. Running pytest..."
if command -v pytest &> /dev/null; then
    if pytest tests/ -q --tb=no 2>/dev/null; then
        check_pass "All pytest tests passed"
    else
        check_fail "Some pytest tests failed"
    fi
else
    check_warn "pytest not found - skipping tests"
fi

# -----------------------------------------------------------------------------
# 7. Rust Build Check
# -----------------------------------------------------------------------------
echo ""
echo "7. Checking Rust build..."
if command -v maturin &> /dev/null; then
    if maturin build --release -m rust/Cargo.toml -q 2>/dev/null; then
        check_pass "Rust builds successfully with maturin"
    else
        check_fail "Rust build failed"
    fi
else
    check_warn "maturin not found - skipping Rust build check"
fi

# -----------------------------------------------------------------------------
# 8. CLI Entry Points Check
# -----------------------------------------------------------------------------
echo ""
echo "8. Checking CLI entry points..."
CLI_COMMANDS=("wasp2-count" "wasp2-map" "wasp2-analyze")
for cmd in "${CLI_COMMANDS[@]}"; do
    if command -v "$cmd" &> /dev/null; then
        if $cmd --help &> /dev/null; then
            check_pass "CLI '$cmd --help' works"
        else
            check_fail "CLI '$cmd --help' failed"
        fi
    else
        check_warn "CLI '$cmd' not installed (run 'pip install -e .' first)"
    fi
done

# -----------------------------------------------------------------------------
# 9. Git Status Check
# -----------------------------------------------------------------------------
echo ""
echo "9. Checking git status..."
UNCOMMITTED=$(git status --porcelain 2>/dev/null | wc -l)
if [ "$UNCOMMITTED" -eq 0 ]; then
    check_pass "Working directory is clean"
else
    check_warn "Found $UNCOMMITTED uncommitted changes"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "VALIDATION SUMMARY"
echo "=============================================="
echo -e "${GREEN}Passed${NC}: $PASSED"
echo -e "${RED}Failed${NC}: $FAILED"
echo ""

if [ "$FAILED" -eq 0 ]; then
    echo -e "${GREEN}All checks passed! Ready to push to upstream.${NC}"
    exit 0
else
    echo -e "${RED}Some checks failed. Please fix before pushing.${NC}"
    exit 1
fi
