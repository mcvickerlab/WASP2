#!/bin/bash
set -e

echo "🚀 Setting up WASP2 development environment..."

# Update conda
echo "📦 Updating conda..."
conda update -n base -c defaults conda -y

# Create environment from environment.yml
echo "🔧 Creating WASP2 conda environment..."
conda env create -f environment.yml -y || conda env update -f environment.yml --prune -y

# Activate environment
echo "✅ Activating WASP2 environment..."
source /opt/conda/etc/profile.d/conda.sh
conda activate WASP2

# Install package in editable mode
echo "📦 Installing WASP2 package in editable mode..."
pip install -e .

# Install development dependencies
echo "🛠️  Installing development dependencies..."
pip install -e ".[dev]" || true

# Set up pre-commit hooks
echo "🪝 Setting up pre-commit hooks..."
pre-commit install || echo "⚠️  Pre-commit not installed, skipping hooks"

# Verify installation
echo "✅ Verifying installation..."
python --version
conda --version
which python

echo ""
echo "✨ WASP2 development environment ready!"
echo ""
echo "📝 Quick start:"
echo "  - Run tests: pytest tests/"
echo "  - Type check: mypy src/"
echo "  - Format code: black src/"
echo "  - Run pipeline: wasp2-count --help"
echo ""
