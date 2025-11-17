# CI/CD Setup for WASP2

**Status**: ✅ Active
**Last Updated**: 2025-11-17
**Based on**: [GenVarLoader](https://github.com/mcvickerlab/GenVarLoader) best practices

---

## Overview

WASP2 uses GitHub Actions for continuous integration and pre-commit hooks for local validation. This ensures:

- ✅ Type hints are validated on every push
- ✅ Regression tests run automatically
- ✅ Full pipeline validation on every PR
- ✅ Code quality checks before commits
- ✅ Multi-Python version testing (3.10, 3.11)

---

## GitHub Actions Workflow

### Triggers
The test workflow (`.github/workflows/test.yml`) runs on:
- **Push** to `main` or any `claude/**` branches
- **Pull requests** targeting `main`
- **Manual trigger** via GitHub UI (workflow_dispatch)

### What Gets Tested

1. **Type Checking** (mypy)
   ```bash
   mypy src/counting/ --ignore-missing-imports
   mypy src/mapping/ --ignore-missing-imports
   mypy src/analysis/ --ignore-missing-imports
   ```

2. **Regression Tests** (pytest)
   ```bash
   pytest tests/regression/ -v --tb=short
   ```
   - Memory regression tests
   - Performance regression tests
   - Output validation tests
   - Full pipeline integration test

3. **Full Pipeline Validation**
   ```bash
   bash scripts/run_full_pipeline_baseline.sh
   ```
   - Counting: 111,454 SNPs
   - Analysis: 43 genomic regions
   - Beta-binomial optimization

4. **Code Coverage**
   ```bash
   pytest --cov=src --cov-report=xml
   ```

### Test Matrix
Tests run on:
- **Python 3.10** (Ubuntu latest)
- **Python 3.11** (Ubuntu latest)

### System Dependencies
Automatically installed:
- `bcftools` - Variant calling
- `bedtools` - Genomic intervals
- `samtools` - BAM file processing
- `time` - Memory/performance profiling

### Python Dependencies
Automatically installed:
- Core: `numpy`, `pandas`, `polars`, `scipy`
- Genomics: `pysam`, `pybedtools`, `anndata`, `scanpy`
- CLI: `typer`, `rich`
- Testing: `pytest`, `pytest-cov`, `mypy`

---

## Pre-commit Hooks

### Installation

```bash
# Install pre-commit
pip install pre-commit

# Install hooks for this repo
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

### What Gets Checked

1. **File Quality**
   - Remove trailing whitespace
   - Fix end-of-file newlines
   - Check YAML syntax
   - Prevent large file commits (>5MB)
   - Detect merge conflicts

2. **Code Formatting** (Black)
   - Line length: 100 characters
   - Consistent Python style

3. **Code Linting** (Flake8)
   - Max line length: 100
   - Ignore: E203, W503 (Black compatibility)

4. **Type Checking** (mypy)
   - Check all `src/` files
   - Ignore missing imports
   - Ensure type hints are valid

5. **Quick Tests** (pytest)
   - Run fast regression tests only (`-m "not slow"`)
   - Catch obvious breakage before commit

### Bypassing Hooks

**Not recommended**, but if absolutely necessary:
```bash
git commit --no-verify -m "message"
```

---

## Local Testing

### Quick Check (before commit)
```bash
# Type checking
mypy src/counting/ src/mapping/ src/analysis/

# Fast tests only
pytest tests/regression/ -v -m "not slow"
```

### Full Validation (before PR)
```bash
# All tests
pytest tests/regression/ -v

# Full pipeline
bash scripts/run_full_pipeline_baseline.sh

# Coverage report
pytest tests/regression/ --cov=src --cov-report=html
open htmlcov/index.html
```

---

## CI/CD Best Practices (from GenVarLoader)

### ✅ What We Implemented

1. **Automated Testing**
   - Multi-Python version matrix
   - System dependencies auto-installed
   - Full pipeline validation

2. **Pre-commit Hooks**
   - Mandatory code quality checks
   - Type validation
   - Quick test runs

3. **PR Requirements**
   - All tests must pass
   - Type checking must pass
   - No option to merge on failure

4. **Test Coverage Tracking**
   - Coverage reports generated
   - Uploaded as artifacts
   - Track coverage trends

5. **Package Validation** ✅ NEW!
   - Test pip installation
   - Verify CLI commands work
   - Build distribution packages
   - Validate with twine

6. **Documentation Build** ✅ NEW!
   - Build Sphinx documentation
   - Check for warnings/errors
   - Validate API autodoc
   - ReadTheDocs configuration

### 📋 Future Enhancements

1. **Automated Versioning**
   - Auto-bump versions on release
   - Semantic versioning

2. **Release Automation**
   - Auto-publish to PyPI
   - GitHub releases with changelogs

4. **Performance Benchmarking**
   - Track performance over time
   - Alert on regressions >30%

---

## Workflow Details

### On Every Push
```
1. Checkout code
2. Set up Python (3.10 & 3.11)
3. Install system deps (bcftools, bedtools, samtools)
4. Install Python deps (including Sphinx, build tools)
5. Verify installations
6. Run mypy type checking ← Validates all type hints!
7. Run regression tests ← Uses our test suite!
8. Run full pipeline ← Real genomic data!
9. Generate coverage report
10. Upload artifacts
11. Test pip installation ← NEW!
12. Build distribution packages ← NEW!
13. Build Sphinx documentation ← NEW!
14. Check docs for warnings ← NEW!
```

### On Every Commit (Pre-commit)
```
1. Format code (Black)
2. Lint code (Flake8)
3. Type check (mypy)
4. Run quick tests
5. Fix file issues
```

---

## Monitoring CI Status

### GitHub UI
- View workflow runs: `Actions` tab in GitHub
- See test results per Python version
- Download coverage reports

### Status Badges
Add to README.md:
```markdown
![Tests](https://github.com/Jaureguy760/WASP2-exp/workflows/WASP2%20Tests/badge.svg)
```

### Email Notifications
GitHub automatically emails on:
- Workflow failures
- First success after failure
- Security vulnerabilities

---

## Troubleshooting

### Pre-commit Hook Failures

**Issue**: Black formatting fails
```bash
# Fix automatically
black src/ tests/ --line-length=100
```

**Issue**: mypy type errors
```bash
# Run mypy locally to see errors
mypy src/analysis/file.py --ignore-missing-imports
```

**Issue**: Quick tests fail
```bash
# Run tests locally with more detail
pytest tests/regression/ -v --tb=long
```

### GitHub Actions Failures

**Issue**: System dependency not found
- Check `.github/workflows/test.yml`
- Ensure `apt-get install` includes the tool

**Issue**: Python import error
- Check pip install section in workflow
- Add missing package

**Issue**: Test timeout
- Increase timeout in pytest
- Mark slow tests with `@pytest.mark.slow`

---

## Files in This Setup

```
.github/
├── workflows/
│   └── test.yml           # Main CI/CD workflow
└── CICD.md                # This documentation

.pre-commit-config.yaml    # Pre-commit hooks config
pytest.ini                 # Pytest configuration
```

---

## Comparison with GenVarLoader

| Feature | GenVarLoader | WASP2 | Status |
|---------|--------------|-------|--------|
| GitHub Actions | ✅ | ✅ | Implemented |
| Pre-commit hooks | ✅ | ✅ | Implemented |
| Multi-Python testing | ✅ (3.10, 3.11, 3.12) | ✅ (3.10, 3.11) | Implemented |
| Type checking in CI | ❓ | ✅ | **Better!** |
| Full pipeline in CI | ❓ | ✅ | **Better!** |
| Automated versioning | ✅ | ❌ | Future work |
| ReadTheDocs | ✅ | ❌ | Future work |

---

## For Contributors

### Before Your First Commit
```bash
# 1. Install pre-commit
pip install pre-commit

# 2. Install hooks
pre-commit install

# 3. Run on existing files
pre-commit run --all-files
```

### Before Opening a PR
```bash
# 1. Ensure type hints are valid
mypy src/counting/ src/mapping/ src/analysis/

# 2. Run full test suite
pytest tests/regression/ -v

# 3. Run full pipeline
bash scripts/run_full_pipeline_baseline.sh

# 4. All green? Open PR!
```

### After PR is Opened
- GitHub Actions will run automatically
- Check the `Actions` tab for results
- Fix any failures and push updates
- PR cannot be merged until all checks pass ✅

---

## Summary

WASP2 now has **production-grade CI/CD** based on GenVarLoader best practices:

✅ **Automated testing** on every push
✅ **Type checking** validates all 24 files
✅ **Pre-commit hooks** catch issues early
✅ **Full pipeline validation** with real genomic data
✅ **Multi-Python testing** ensures compatibility
✅ **Coverage tracking** monitors test quality

**Result**: High-confidence deployments and easy collaboration! 🚀
