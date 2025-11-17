# WASP2 Development Setup

This guide covers setting up your development environment for WASP2.

## Quick Start

### Option 1: Using Conda/Mamba (Recommended)

```bash
# Create environment from environment.yml
conda env create -f environment.yml

# Activate environment
conda activate WASP2

# Verify installation
pytest tests/
```

### Option 2: Using pip + system packages

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install bedtools bcftools samtools

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt

# Verify installation
pytest tests/
```

## Running Tests

### Full Test Suite
```bash
python -m pytest tests/ -v
```

### Regression Tests Only
```bash
python -m pytest tests/regression/ -v
```

### With Coverage
```bash
pytest tests/ --cov=src --cov-report=html
```

## Type Checking

```bash
# Check all modules
mypy src/

# Check specific module
mypy src/counting/
mypy src/mapping/
mypy src/analysis/
```

## Project Structure

```
WASP2-exp/
├── src/
│   ├── counting/      # Allele counting module (type hints ✓)
│   ├── mapping/       # WASP mapping module (type hints ✓)
│   └── analysis/      # Statistical analysis module
├── tests/
│   └── regression/    # Regression test suite
├── test_data/         # Test data and baselines
├── mypy.ini           # Type checking configuration
└── environment.yml    # Conda dependencies
```

## Development Workflow

1. **Make changes** to source code
2. **Run type checker**: `mypy src/`
3. **Run tests**: `pytest tests/`
4. **Commit changes**: Following conventional commits

## Type Hints Status

| Module   | Files | Lines | Status |
|----------|-------|-------|--------|
| Counting | 7     | 1,424 | ✅ Complete |
| Mapping  | 7     | 1,569 | ✅ Complete |
| Analysis | 10    | 2,507 | 📋 Planned |

## Common Issues

### Tests fail with "ModuleNotFoundError"
- **Solution**: Ensure you're using `python -m pytest` instead of just `pytest`
- This ensures tests run with the correct Python interpreter

### mypy errors about missing stubs
- **Expected**: Third-party libraries (pysam, polars, scipy) don't have type stubs
- These are configured to be ignored in `mypy.ini`

### Type hints cause runtime errors
- **Not expected**: Type hints are annotations only and don't affect runtime
- If this occurs, it's a separate bug - file an issue!
