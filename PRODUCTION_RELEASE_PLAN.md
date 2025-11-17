# WASP2 Production Release Plan

**Goal**: Package WASP2 for PyPI and create professional documentation, then merge to main as v1.0.0

**Status**: 📋 Planning Phase
**Started**: 2025-11-17
**Target Release**: v1.0.0

---

## Phase 1: PyPI Package Setup (Option 3)

### Objectives
- Make WASP2 `pip install`-able
- Professional package metadata
- CLI commands accessible system-wide
- Dependencies properly declared

### Tasks

#### 1.1 Create `pyproject.toml` (Modern Python Packaging)
**Why `pyproject.toml` not `setup.py`?**
- Modern Python standard (PEP 518, 621)
- Cleaner than setup.py
- Better tool integration
- Required by Poetry, Hatch, etc.

**Contents**:
```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "wasp2"
version = "1.0.0"
description = "Allele-specific analysis of next-gen sequencing data"
authors = [{name = "McVicker Lab"}]
license = {text = "MIT"}
readme = "README.md"
requires-python = ">=3.10"
classifiers = [...]
dependencies = [
    "numpy>=1.21.0",
    "pandas>=2.0.0",
    "polars>=0.19.0",
    "scipy>=1.10.0",
    "pysam>=0.21.0",
    "pybedtools>=0.9.0",
    "anndata>=0.8.0",
    "typer>=0.9.0",
    "rich>=13.0.0"
]

[project.optional-dependencies]
dev = ["pytest>=7.0", "mypy>=1.0", "black>=23.0", "pre-commit>=3.0"]

[project.scripts]
wasp2-count = "counting.__main__:app"
wasp2-map = "mapping.__main__:app"
wasp2-analyze = "analysis.__main__:app"

[project.urls]
Homepage = "https://github.com/Jaureguy760/WASP2-exp"
Documentation = "https://wasp2.readthedocs.io"
Repository = "https://github.com/Jaureguy760/WASP2-exp"
```

#### 1.2 Fix Source Layout
**Current structure**:
```
WASP2-exp/
├── src/
│   ├── counting/
│   ├── mapping/
│   └── analysis/
```

**Issue**: Need to make `src/` a proper package

**Options**:

**Option A: Keep flat structure** (simpler)
```python
# Users import like:
from counting import count_alleles
from mapping import remap_utils
from analysis import as_analysis
```

**Option B: Unified package** (cleaner)
```
src/
└── wasp2/
    ├── __init__.py
    ├── counting/
    ├── mapping/
    └── analysis/
```

```python
# Users import like:
from wasp2.counting import count_alleles
from wasp2.mapping import remap_utils
from wasp2.analysis import as_analysis
```

**Recommendation**: Option B (unified package)
- Cleaner imports
- Avoids namespace conflicts
- Standard Python practice

**Migration needed**:
1. Move src/* to src/wasp2/*
2. Add src/wasp2/__init__.py
3. Update all imports in code
4. Update CLI entry points

#### 1.3 Create MANIFEST.in
Include non-Python files in package:
```
include README.md
include LICENSE
include CITATION.cff
recursive-include test_data *.bam *.vcf *.bed
recursive-include scripts *.sh
```

#### 1.4 Test Local Installation
```bash
# Editable install (for development)
pip install -e .

# Test CLI commands work
wasp2-count --help
wasp2-map --help
wasp2-analyze --help

# Test imports
python -c "from wasp2.counting import count_alleles"
python -c "from wasp2.mapping import remap_utils"
python -c "from wasp2.analysis import as_analysis"

# Run tests with installed package
pytest tests/
```

#### 1.5 Build Package Locally
```bash
# Install build tools
pip install build twine

# Build distribution
python -m build

# Verify built package
ls dist/
# Should see:
# wasp2-1.0.0.tar.gz
# wasp2-1.0.0-py3-none-any.whl

# Check package
twine check dist/*
```

---

## Phase 2: Documentation Site (Option 4)

### Objectives
- Professional docs on ReadTheDocs
- API documentation auto-generated from type hints
- Installation guide
- Tutorials and examples

### Tasks

#### 2.1 Set Up Sphinx
```bash
# Install Sphinx
pip install sphinx sphinx-rtd-theme sphinx-autodoc-typehints

# Initialize docs
mkdir docs
cd docs
sphinx-quickstart

# Choose:
# - Separate source/build: yes
# - autodoc: yes
# - intersphinx: yes
```

#### 2.2 Configure Sphinx (`docs/source/conf.py`)
```python
import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

project = 'WASP2'
author = 'McVicker Lab'
version = '1.0.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # Google/NumPy docstring support
    'sphinx.ext.viewcode',   # Add source code links
    'sphinx.ext.intersphinx', # Link to other docs
    'sphinx_autodoc_typehints', # Use our type hints!
]

html_theme = 'sphinx_rtd_theme'
```

#### 2.3 Documentation Structure
```
docs/
├── source/
│   ├── index.rst              # Landing page
│   ├── installation.rst       # pip install wasp2
│   ├── quickstart.rst         # 5-minute tutorial
│   ├── user_guide/
│   │   ├── counting.rst       # Count alleles
│   │   ├── mapping.rst        # WASP remapping
│   │   └── analysis.rst       # Statistical analysis
│   ├── api/
│   │   ├── counting.rst       # Auto-generated API
│   │   ├── mapping.rst        # Auto-generated API
│   │   └── analysis.rst       # Auto-generated API
│   ├── development/
│   │   ├── contributing.rst   # How to contribute
│   │   ├── testing.rst        # Running tests
│   │   └── cicd.rst           # CI/CD pipeline
│   └── changelog.rst          # Version history
├── Makefile
└── make.bat
```

#### 2.4 Create Key Documentation Pages

**index.rst** (Landing page):
```rst
WASP2: Allele-Specific Analysis
================================

.. image:: https://img.shields.io/pypi/v/wasp2
.. image:: https://github.com/Jaureguy760/WASP2-exp/workflows/WASP2%20Tests/badge.svg

WASP2 is a suite of tools for unbiased allele-specific analysis of
next-generation sequencing data.

Installation
------------

.. code-block:: bash

   pip install wasp2

Quick Example
-------------

.. code-block:: python

   from wasp2.counting import count_alleles

   # Count alleles from BAM file
   counts = count_alleles("sample.bam", "variants.vcf")

Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   user_guide/index
   api/index
   development/index
   changelog
```

**installation.rst**:
```rst
Installation
============

Requirements
------------

System Dependencies
~~~~~~~~~~~~~~~~~~~

WASP2 requires the following tools:

- bcftools >= 1.10
- bedtools >= 2.29
- samtools >= 1.10

On Ubuntu/Debian:

.. code-block:: bash

   sudo apt-get install bcftools bedtools samtools

Python Requirements
~~~~~~~~~~~~~~~~~~~

Python 3.10 or higher is required.

Install via pip
---------------

.. code-block:: bash

   pip install wasp2

Development Installation
------------------------

.. code-block:: bash

   git clone https://github.com/Jaureguy760/WASP2-exp
   cd WASP2-exp
   pip install -e ".[dev]"
```

**quickstart.rst**:
```rst
Quick Start
===========

This 5-minute tutorial shows basic WASP2 usage.

Count Alleles
-------------

.. code-block:: bash

   wasp2-count count-alleles \\
     --bam sample.bam \\
     --vcf variants.vcf \\
     --output counts.tsv

WASP Mapping
------------

.. code-block:: bash

   wasp2-map find-intersecting-snps \\
     --bam input.bam \\
     --vcf variants.vcf \\
     --output intersecting.bam

Analyze Imbalance
-----------------

.. code-block:: bash

   wasp2-analyze find-imbalance \\
     --count-file counts.tsv \\
     --output results.tsv
```

#### 2.5 Auto-Generate API Docs from Type Hints

**api/counting.rst**:
```rst
Counting Module API
===================

.. automodule:: wasp2.counting.count_alleles
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: wasp2.counting.filter_variant_data
   :members:
   :undoc-members:
   :show-inheritance:
```

**KEY**: Because we have comprehensive type hints, Sphinx will auto-generate:
- Function signatures with types
- Parameter descriptions
- Return types
- Example usage

#### 2.6 Configure ReadTheDocs

Create `.readthedocs.yaml`:
```yaml
version: 2

build:
  os: ubuntu-22.04
  tools:
    python: "3.11"

sphinx:
  configuration: docs/source/conf.py

python:
  install:
    - method: pip
      path: .
      extra_requirements:
        - dev
```

#### 2.7 Test Docs Build
```bash
cd docs
make html

# Open in browser
open build/html/index.html

# Check for warnings
make clean html 2>&1 | grep -i warning
```

---

## Phase 3: Update CI/CD

### Add to `.github/workflows/test.yml`

```yaml
- name: Test package installation
  run: |
    pip install -e .
    wasp2-count --version
    wasp2-map --version
    wasp2-analyze --version

- name: Build package
  run: |
    pip install build
    python -m build
    twine check dist/*

- name: Build documentation
  run: |
    cd docs
    make html

- name: Check docs for warnings
  run: |
    cd docs
    make clean html 2>&1 | tee build.log
    if grep -i warning build.log; then
      echo "Documentation has warnings!"
      exit 1
    fi
```

---

## Phase 4: Create Comprehensive PR

### PR Title
```
Production Release: Type Hints (100%) + PyPI Package + Documentation
```

### PR Description

```markdown
## 🎉 WASP2 v1.0.0 Production Release

This PR transforms WASP2 into a production-ready, pip-installable package with comprehensive documentation.

### Type Hints Coverage (100%) ✅

**Statistics:**
- 24 files fully typed
- 5,500 lines of code with type annotations
- 0 mypy errors

**Modules:**
- ✅ Counting: 7 files (1,424 lines)
- ✅ Mapping: 7 files (1,569 lines)
- ✅ Analysis: 10 files (2,507 lines)

**Validation:**
- All type checks pass
- Full pipeline tested with real genomic data
- Memory & performance regression tests passing

### PyPI Package Setup ✅

**Installation:**
```bash
pip install wasp2
```

**Features:**
- Modern pyproject.toml configuration
- CLI commands: `wasp2-count`, `wasp2-map`, `wasp2-analyze`
- Proper dependency management
- Development extras: `pip install wasp2[dev]`

### Documentation Site ✅

**ReadTheDocs:** https://wasp2.readthedocs.io

**Contents:**
- Installation guide
- Quick start tutorial
- User guides for each module
- **API documentation auto-generated from type hints**
- Development guide
- CI/CD documentation

### CI/CD Pipeline ✅

**GitHub Actions:**
- Runs on every push/PR
- Multi-Python testing (3.10, 3.11)
- Type checking with mypy
- Full regression test suite
- Package build validation
- Documentation build validation

**Pre-commit hooks:**
- Code formatting (Black)
- Linting (Flake8)
- Type checking (mypy)
- Quick tests

### Commits
1. Bug fixes and initial typing
2. TH-3: Mapping module type hints
3. TH-2 Wave 1-5: Analysis module type hints
4. CI/CD setup
5. Full pipeline validation
6. PyPI package configuration
7. Sphinx documentation setup

### Testing
All tests passing:
- ✅ 24 files pass mypy type checking
- ✅ Regression tests (memory, performance)
- ✅ Full pipeline with real genomic data
- ✅ Package installation
- ✅ Documentation builds without warnings

### Breaking Changes
None - this is a new major version with added features.

### Migration Guide
For existing users:
```bash
# Old way (from source)
python -m counting.count_alleles ...

# New way (pip installed)
wasp2-count count-alleles ...

# Or programmatically
from wasp2.counting import count_alleles
```

### Next Steps After Merge
1. Tag v1.0.0 release
2. Publish to PyPI
3. Announce release
4. Update README with pip install instructions
```

---

## Phase 5: Release & Publish

### After PR is Merged to Main

```bash
# 1. Update local main
git checkout main
git pull origin main

# 2. Tag release
git tag -a v1.0.0 -m "WASP2 v1.0.0 - Production Release

- Complete type hint coverage (24 files)
- PyPI package available
- Professional documentation on ReadTheDocs
- CI/CD pipeline with GitHub Actions
- Pre-commit hooks
"

# 3. Push tag
git push origin v1.0.0

# 4. Build package
python -m build

# 5. Publish to PyPI (REQUIRES API TOKEN)
twine upload dist/*

# 6. Create GitHub release
gh release create v1.0.0 \
  --title "WASP2 v1.0.0 - Production Release" \
  --notes "See CHANGELOG.md for details" \
  dist/*
```

---

## Timeline Estimate

| Phase | Tasks | Time Estimate |
|-------|-------|---------------|
| **1. PyPI Package** | pyproject.toml, imports, testing | 2-3 hours |
| **2. Documentation** | Sphinx setup, write docs, API | 3-4 hours |
| **3. Update CI/CD** | Add package/docs checks | 1 hour |
| **4. Create PR** | Write comprehensive PR | 1 hour |
| **5. Release** | Tag, publish, announce | 1 hour |
| **TOTAL** | | **8-10 hours** |

---

## Success Criteria

✅ **Package Works:**
- `pip install wasp2` succeeds
- CLI commands work: `wasp2-count`, `wasp2-map`, `wasp2-analyze`
- Imports work: `from wasp2.counting import count_alleles`
- All dependencies resolved

✅ **Documentation Live:**
- ReadTheDocs builds successfully
- All pages render correctly
- API docs show type hints
- Examples work

✅ **CI/CD Passes:**
- All tests pass
- Package builds successfully
- Docs build without warnings
- Type checking passes

✅ **Release Published:**
- v1.0.0 tag created
- Published on PyPI
- GitHub release created
- README updated

---

## Risk Assessment

### Medium Risk: Import Refactoring

**Issue**: Moving to `wasp2.*` package structure requires updating all internal imports

**Mitigation**:
1. Use find/replace carefully
2. Test each module after changes
3. Run full test suite
4. mypy will catch import errors

### Low Risk: PyPI Name Availability

**Issue**: `wasp2` might be taken on PyPI

**Check**:
```bash
pip search wasp2  # or check https://pypi.org/project/wasp2/
```

**Mitigation**:
- Alternative names: `wasp-bio`, `wasp2-analysis`, `wasp2-mcvickerlab`

### Low Risk: Documentation Build Issues

**Issue**: Sphinx might have warnings or errors

**Mitigation**:
- Test locally first
- Fix all warnings before pushing
- CI will catch build failures

---

## Rollback Plan

If anything goes wrong:

1. **PR not yet merged**: Close PR, fix issues, reopen
2. **Merged but not released**: Revert merge commit
3. **Released to PyPI**: Cannot delete, must release v1.0.1 with fixes

**Prevention**: Thorough testing before each phase!

---

## Post-Release Tasks

1. **Update README.md**:
   - Add pip install instructions
   - Add PyPI badge
   - Add ReadTheDocs badge

2. **Announce**:
   - Post to lab Slack/email
   - Tweet (if applicable)
   - Update any related papers

3. **Monitor**:
   - Watch for PyPI download stats
   - Check ReadTheDocs analytics
   - Monitor GitHub issues

4. **Plan v1.1.0**:
   - Gather feedback
   - Plan new features
   - Continue improvements

---

## Questions to Resolve

1. **Package name**: `wasp2` available on PyPI?
2. **License**: MIT? GPL? Apache?
3. **Authors**: Who to credit in pyproject.toml?
4. **Repository**: Keep Jaureguy760/WASP2-exp or move to mcvickerlab org?
5. **PyPI account**: Who has access to publish?

---

## Ready to Start!

Phase 1 begins now: Setting up PyPI package structure.

**First task**: Create `pyproject.toml` with all metadata.
