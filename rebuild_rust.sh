#!/bin/bash
# Rebuild Rust extension with indel support
# This script rebuilds the Rust filter with same-locus slop parameter

set -e

echo "🔧 Rebuilding WASP2 Rust extension with indel support..."

# Set LIBCLANG_PATH
export LIBCLANG_PATH=/iblm/netapp/home/jjaureguy/mambaforge/lib/python3.10/site-packages/clang/native
export LD_LIBRARY_PATH=/iblm/netapp/home/jjaureguy/mambaforge/lib:$LD_LIBRARY_PATH

# Navigate to rust directory
cd rust

# Clean previous build
echo "📦 Cleaning previous build..."
cargo clean

# Build with maturin
echo "🦀 Building Rust extension..."
maturin develop --release

echo "✅ Rust extension rebuilt successfully!"
echo ""
echo "Test it with:"
echo "  python -c \"from wasp2_rust import filter_bam_wasp; import inspect; print(inspect.signature(filter_bam_wasp))\""
echo ""
echo "Expected output should include: same_locus_slop=0"
