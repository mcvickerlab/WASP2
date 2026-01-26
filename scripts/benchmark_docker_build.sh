#!/bin/bash
# WASP2 Docker Build Benchmark Script
# Compares original vs optimized Dockerfile builds
# Usage: ./scripts/benchmark_docker_build.sh

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

echo "========================================"
echo "WASP2 Docker Build Benchmark"
echo "========================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to format time
format_time() {
    local seconds=$1
    printf "%dm %ds" $((seconds/60)) $((seconds%60))
}

# Function to format size
format_size() {
    local size=$1
    if [[ $size =~ ^[0-9]+$ ]]; then
        printf "%.1fMB" $(echo "scale=1; $size/1048576" | bc)
    else
        echo "$size"
    fi
}

# Ensure BuildKit is enabled
export DOCKER_BUILDKIT=1

echo "Build Configuration:"
echo "  - BuildKit: enabled"
echo "  - Project: $PROJECT_DIR"
echo ""

# Benchmark function
benchmark_build() {
    local name=$1
    local dockerfile=$2
    local tag=$3

    echo "----------------------------------------"
    echo "Building: $name"
    echo "  Dockerfile: $dockerfile"
    echo "  Tag: $tag"
    echo "----------------------------------------"

    # Clear build cache for fair comparison
    docker builder prune -f --filter type=exec.cachemount 2>/dev/null || true

    # Build with timing
    local start_time=$(date +%s)

    if docker buildx build \
        -f "$dockerfile" \
        -t "$tag" \
        --progress=plain \
        --no-cache \
        . 2>&1 | tee /tmp/build_${name}.log; then

        local end_time=$(date +%s)
        local duration=$((end_time - start_time))

        # Get image size
        local size=$(docker images "$tag" --format "{{.Size}}")

        echo ""
        echo -e "${GREEN}SUCCESS:${NC} $name"
        echo "  Build time: $(format_time $duration)"
        echo "  Image size: $size"

        # Store results
        echo "$name,$duration,$size" >> /tmp/benchmark_results.csv
    else
        echo -e "${RED}FAILED:${NC} $name"
        echo "$name,FAILED,FAILED" >> /tmp/benchmark_results.csv
    fi
    echo ""
}

# Benchmark cached rebuild
benchmark_cached_build() {
    local name=$1
    local dockerfile=$2
    local tag=$3

    echo "----------------------------------------"
    echo "Cached Rebuild: $name"
    echo "----------------------------------------"

    # Make a small change to trigger rebuild
    touch src/counting/__init__.py

    local start_time=$(date +%s)

    if docker buildx build \
        -f "$dockerfile" \
        -t "$tag" \
        --progress=plain \
        . 2>&1 | tee /tmp/cached_${name}.log; then

        local end_time=$(date +%s)
        local duration=$((end_time - start_time))

        echo -e "${GREEN}Cached rebuild time:${NC} $(format_time $duration)"
        echo "$name-cached,$duration,-" >> /tmp/benchmark_results.csv
    fi
    echo ""
}

# Initialize results file
echo "build,duration_seconds,size" > /tmp/benchmark_results.csv

# Check if both Dockerfiles exist
if [[ ! -f "Dockerfile" ]]; then
    echo -e "${RED}ERROR:${NC} Original Dockerfile not found"
    exit 1
fi

if [[ ! -f "Dockerfile.optimized" ]]; then
    echo -e "${YELLOW}WARNING:${NC} Dockerfile.optimized not found, skipping optimized build"
    SKIP_OPTIMIZED=1
fi

# Run benchmarks
echo ""
echo "========================================"
echo "Phase 1: Fresh Builds (no cache)"
echo "========================================"

benchmark_build "original" "Dockerfile" "wasp2:original"

if [[ -z "$SKIP_OPTIMIZED" ]]; then
    benchmark_build "optimized" "Dockerfile.optimized" "wasp2:optimized"
fi

echo ""
echo "========================================"
echo "Phase 2: Cached Rebuilds"
echo "========================================"

benchmark_cached_build "original" "Dockerfile" "wasp2:original"

if [[ -z "$SKIP_OPTIMIZED" ]]; then
    benchmark_cached_build "optimized" "Dockerfile.optimized" "wasp2:optimized"
fi

echo ""
echo "========================================"
echo "Results Summary"
echo "========================================"
echo ""
cat /tmp/benchmark_results.csv | column -t -s ','
echo ""

# Security scan comparison
echo "========================================"
echo "Security Scan (Trivy)"
echo "========================================"

if command -v trivy &> /dev/null; then
    echo ""
    echo "Scanning wasp2:original..."
    trivy image --severity HIGH,CRITICAL --quiet wasp2:original 2>/dev/null || echo "Scan failed or no issues"

    if [[ -z "$SKIP_OPTIMIZED" ]]; then
        echo ""
        echo "Scanning wasp2:optimized..."
        trivy image --severity HIGH,CRITICAL --quiet wasp2:optimized 2>/dev/null || echo "Scan failed or no issues"
    fi
else
    echo "Trivy not installed. Install with: brew install trivy"
fi

echo ""
echo "========================================"
echo "Benchmark Complete"
echo "========================================"
echo ""
echo "Full build logs available at:"
echo "  /tmp/build_original.log"
if [[ -z "$SKIP_OPTIMIZED" ]]; then
    echo "  /tmp/build_optimized.log"
fi
