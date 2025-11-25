#!/bin/bash
# =============================================================================
# WASP2 Development Environment Setup Script
# =============================================================================
#
# This script sets up the complete development environment for WASP2,
# including all dependencies for multi-format variant support (VCF + PLINK2)
#
# Usage:
#   ./scripts/setup_dev_env.sh [--full]
#
# Options:
#   --full    Install all optional dependencies (docs, dev tools)
#
# Requirements:
#   - Conda or Mamba installed
#   - Rust toolchain installed (for Rust extension)
#
# =============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Parse arguments
FULL_INSTALL=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --full)
            FULL_INSTALL=true
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  WASP2 Development Environment Setup  ${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "Project root: ${GREEN}$PROJECT_ROOT${NC}"
echo ""

# -----------------------------------------------------------------------------
# Step 1: Check prerequisites
# -----------------------------------------------------------------------------
echo -e "${YELLOW}Step 1: Checking prerequisites...${NC}"

# Check for conda/mamba
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo -e "  ✓ Found mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    echo -e "  ✓ Found conda"
else
    echo -e "${RED}  ✗ Neither conda nor mamba found. Please install one first.${NC}"
    exit 1
fi

# Check for Rust
if command -v rustc &> /dev/null; then
    RUST_VERSION=$(rustc --version)
    echo -e "  ✓ Found Rust: $RUST_VERSION"
else
    echo -e "${RED}  ✗ Rust not found. Please install Rust first:${NC}"
    echo -e "    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python --version 2>&1 | cut -d' ' -f2)
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d'.' -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d'.' -f2)
if [[ $PYTHON_MAJOR -eq 3 && $PYTHON_MINOR -ge 10 ]]; then
    echo -e "  ✓ Python version: $PYTHON_VERSION"
else
    echo -e "${RED}  ✗ Python 3.10+ required, found $PYTHON_VERSION${NC}"
    exit 1
fi

echo ""

# -----------------------------------------------------------------------------
# Step 2: Install bioinformatics tools via conda
# -----------------------------------------------------------------------------
echo -e "${YELLOW}Step 2: Installing bioinformatics tools...${NC}"

# Check and install plink2
if ! command -v plink2 &> /dev/null; then
    echo -e "  Installing plink2..."
    $CONDA_CMD install -y -c bioconda plink2
    echo -e "  ✓ Installed plink2"
else
    echo -e "  ✓ plink2 already installed"
fi

# Check and install bedtools
if ! command -v bedtools &> /dev/null; then
    echo -e "  Installing bedtools..."
    $CONDA_CMD install -y -c bioconda bedtools
    echo -e "  ✓ Installed bedtools"
else
    echo -e "  ✓ bedtools already installed"
fi

# Check and upgrade htslib/samtools if old
HTSLIB_VERSION=$($CONDA_CMD list htslib 2>/dev/null | grep -v "^#" | awk '{print $2}' || echo "0")
if [[ "$HTSLIB_VERSION" < "1.17" ]]; then
    echo -e "  Upgrading htslib and samtools..."
    $CONDA_CMD install -y -c bioconda 'htslib>=1.17' 'samtools>=1.17'
    echo -e "  ✓ Upgraded htslib/samtools"
else
    echo -e "  ✓ htslib version OK: $HTSLIB_VERSION"
fi

echo ""

# -----------------------------------------------------------------------------
# Step 3: Install Python dependencies
# -----------------------------------------------------------------------------
echo -e "${YELLOW}Step 3: Installing Python dependencies...${NC}"

cd "$PROJECT_ROOT"

# Install core dependencies
echo -e "  Installing core dependencies..."
pip install -q numpy pandas polars scipy pysam pybedtools anndata typer rich

# Install scanpy (can be slow)
if ! python -c "import scanpy" 2>/dev/null; then
    echo -e "  Installing scanpy (this may take a minute)..."
    pip install -q "scanpy>=1.9.0"
    echo -e "  ✓ Installed scanpy"
else
    echo -e "  ✓ scanpy already installed"
fi

# Install pgenlib for PLINK2 support
if ! python -c "import pgenlib" 2>/dev/null; then
    echo -e "  Installing pgenlib..."
    pip install -q Pgenlib
    echo -e "  ✓ Installed pgenlib"
else
    echo -e "  ✓ pgenlib already installed"
fi

# Install dev dependencies
echo -e "  Installing dev dependencies..."
pip install -q pytest pytest-cov maturin build twine

if [ "$FULL_INSTALL" = true ]; then
    echo -e "  Installing full dev tools..."
    pip install -q black flake8 mypy pre-commit

    echo -e "  Installing documentation tools..."
    pip install -q "sphinx>=5.0" "pydata-sphinx-theme>=0.14" "sphinx-autodoc-typehints>=1.0"
fi

echo ""

# -----------------------------------------------------------------------------
# Step 4: Build and install WASP2 in development mode
# -----------------------------------------------------------------------------
echo -e "${YELLOW}Step 4: Building and installing WASP2...${NC}"

cd "$PROJECT_ROOT"

# Build Rust extension and install in development mode
echo -e "  Building Rust extension and installing package..."
if [ -f "rust/Cargo.toml" ]; then
    maturin develop --release -m rust/Cargo.toml 2>&1 | tail -5
    echo -e "  ✓ Built Rust extension"
else
    echo -e "${YELLOW}  ⚠ No Rust extension found, installing Python-only...${NC}"
    pip install -e .
fi

echo ""

# -----------------------------------------------------------------------------
# Step 5: Verify installation
# -----------------------------------------------------------------------------
echo -e "${YELLOW}Step 5: Verifying installation...${NC}"

# Test core imports
echo -e "  Testing imports..."

IMPORT_ERRORS=0

# Test Rust extension
if python -c "import wasp2_rust" 2>/dev/null; then
    echo -e "  ✓ wasp2_rust (Rust extension)"
else
    echo -e "${RED}  ✗ wasp2_rust failed to import${NC}"
    IMPORT_ERRORS=$((IMPORT_ERRORS + 1))
fi

# Test counting module
if python -c "from counting import __main__" 2>/dev/null; then
    echo -e "  ✓ counting module"
else
    echo -e "${RED}  ✗ counting module failed to import${NC}"
    IMPORT_ERRORS=$((IMPORT_ERRORS + 1))
fi

# Test mapping module
if python -c "from mapping import __main__" 2>/dev/null; then
    echo -e "  ✓ mapping module"
else
    echo -e "${RED}  ✗ mapping module failed to import${NC}"
    IMPORT_ERRORS=$((IMPORT_ERRORS + 1))
fi

# Test analysis module
if python -c "from analysis import __main__" 2>/dev/null; then
    echo -e "  ✓ analysis module"
else
    echo -e "${RED}  ✗ analysis module failed to import${NC}"
    IMPORT_ERRORS=$((IMPORT_ERRORS + 1))
fi

# Test new io module
if python -c "from wasp2.io import VariantSource, VCFSource" 2>/dev/null; then
    echo -e "  ✓ wasp2.io module (VariantSource)"
else
    echo -e "${RED}  ✗ wasp2.io module failed to import${NC}"
    IMPORT_ERRORS=$((IMPORT_ERRORS + 1))
fi

# Test pgenlib
if python -c "import pgenlib" 2>/dev/null; then
    echo -e "  ✓ pgenlib (PLINK2 support)"
else
    echo -e "${YELLOW}  ⚠ pgenlib not available (PGEN support disabled)${NC}"
fi

echo ""

# Test external tools
echo -e "  Testing external tools..."

if command -v bcftools &> /dev/null; then
    echo -e "  ✓ bcftools: $(bcftools --version | head -1)"
else
    echo -e "${RED}  ✗ bcftools not found${NC}"
fi

if command -v bedtools &> /dev/null; then
    echo -e "  ✓ bedtools: $(bedtools --version)"
else
    echo -e "${RED}  ✗ bedtools not found${NC}"
fi

if command -v plink2 &> /dev/null; then
    echo -e "  ✓ plink2: $(plink2 --version | head -1)"
else
    echo -e "${YELLOW}  ⚠ plink2 not found (PGEN conversion disabled)${NC}"
fi

echo ""

# -----------------------------------------------------------------------------
# Step 6: Run quick tests
# -----------------------------------------------------------------------------
echo -e "${YELLOW}Step 6: Running quick tests...${NC}"

cd "$PROJECT_ROOT"
python -m pytest tests/io/test_variant_source.py::TestVariant tests/io/test_variant_source.py::TestGenotype -v --tb=short 2>&1 | tail -15

echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}           Setup Complete!             ${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

if [ $IMPORT_ERRORS -eq 0 ]; then
    echo -e "${GREEN}✓ All modules installed successfully!${NC}"
else
    echo -e "${YELLOW}⚠ $IMPORT_ERRORS module(s) failed to import${NC}"
fi

echo ""
echo -e "You can now use WASP2:"
echo -e "  ${GREEN}wasp2-count --help${NC}    # Allele counting"
echo -e "  ${GREEN}wasp2-map --help${NC}      # WASP mapping"
echo -e "  ${GREEN}wasp2-analyze --help${NC}  # Allelic imbalance analysis"
echo ""
echo -e "Or use the new multi-format API:"
echo -e "  ${GREEN}from wasp2.io import VariantSource${NC}"
echo -e "  ${GREEN}with VariantSource.open('variants.vcf.gz') as src:${NC}"
echo -e "  ${GREEN}    src.to_bed('output.bed')${NC}"
echo ""
echo -e "Run full test suite with:"
echo -e "  ${GREEN}pytest tests/ -v${NC}"
echo ""
