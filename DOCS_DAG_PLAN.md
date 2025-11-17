# Phase 2: Sphinx Documentation DAG & Parallelization

**Objective**: Create professional documentation with API auto-generation in 1.5-2 hours (vs 3-4 sequential)

**Time Savings**: 50% reduction through parallelization

---

## Documentation DAG

```
┌─────────────────────────────────────────────────────────┐
│  Wave 1: Foundation (Sequential - 1 agent)              │
│  Duration: 30 minutes                                   │
├─────────────────────────────────────────────────────────┤
│  - Set up Sphinx directory structure                    │
│  - Create conf.py with extensions                       │
│  - Create basic index.rst                               │
│  - Create Makefile                                      │
└─────────────────┬───────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────┐
│  Wave 2: Content Creation (4 Parallel Agents)           │
│  Duration: ~1 hour (longest agent)                      │
├──────────────┬──────────────┬──────────────┬───────────┤
│  Agent A     │  Agent B     │  Agent C     │  Agent D  │
│  (Basics)    │  (Counting)  │  (Mapping)   │ (Analysis)│
├──────────────┼──────────────┼──────────────┼───────────┤
│ installation │ user_guide/  │ user_guide/  │user_guide/│
│ .rst         │ counting.rst │ mapping.rst  │analysis.  │
│              │              │              │rst        │
│ quickstart   │              │              │           │
│ .rst         │              │              │development│
│              │              │              │.rst       │
└──────┬───────┴──────┬───────┴──────┬───────┴─────┬─────┘
       │              │              │             │
       └──────────────┴──────────────┴─────────────┘
                                 ▼
┌─────────────────────────────────────────────────────────┐
│  Wave 3: API Autodoc (3 Parallel Agents)                │
│  Duration: ~45 minutes (longest agent)                  │
├──────────────────┬──────────────────┬──────────────────┤
│  Agent E         │  Agent F         │  Agent G         │
│  (Counting API)  │  (Mapping API)   │  (Analysis API)  │
├──────────────────┼──────────────────┼──────────────────┤
│ api/counting.rst │ api/mapping.rst  │ api/analysis.rst │
│ (autodoc for     │ (autodoc for     │ (autodoc for     │
│  7 files)        │  7 files)        │  10 files)       │
└──────────┬───────┴──────────┬───────┴──────────┬───────┘
           │                  │                  │
           └──────────────────┴──────────────────┘
                              ▼
┌─────────────────────────────────────────────────────────┐
│  Wave 4: Integration (Sequential - 1 agent)             │
│  Duration: 30 minutes                                   │
├─────────────────────────────────────────────────────────┤
│  - Update index.rst with all toctrees                   │
│  - Create .readthedocs.yaml                             │
│  - Add navigation structure                             │
│  - Build and test locally                               │
│  - Fix any build warnings                               │
└─────────────────────────────────────────────────────────┘
```

---

## Time Analysis

### Sequential Approach (Original)
```
Sphinx setup:        30 min
installation.rst:    20 min
quickstart.rst:      30 min
counting guide:      30 min
mapping guide:       30 min
analysis guide:      40 min
counting API:        25 min
mapping API:         25 min
analysis API:        35 min
integration:         30 min
testing/fixes:       20 min
────────────────────────────
TOTAL: 3h 15min
```

### Parallel Approach (Optimized)
```
Wave 1 (sequential):    30 min
Wave 2 (parallel):      60 min (longest agent)
Wave 3 (parallel):      45 min (longest agent)
Wave 4 (sequential):    30 min
────────────────────────────
TOTAL: 2h 45min - 3h 5min

Savings: 10-30 minutes (but safer & better quality)
```

**Note**: Actually the time is similar, BUT we get:
- Better quality (specialized agents per topic)
- Error isolation (one agent fails ≠ all fail)
- Easier debugging
- Can re-run individual agents if needed

---

## Wave 1: Foundation

**Single Agent Task**

Create basic Sphinx infrastructure:

### Directory Structure
```
docs/
├── source/
│   ├── conf.py           # Sphinx configuration
│   ├── index.rst         # Landing page
│   ├── _static/          # Custom CSS/images
│   ├── _templates/       # Custom templates
│   └── requirements.txt  # Docs dependencies
├── Makefile              # Build commands
└── make.bat              # Windows build
```

### conf.py Configuration
```python
import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

project = 'WASP2'
copyright = '2025, Aaron Ho, Jeff Jaureguy, McVicker Lab'
author = 'Aaron Ho, Jeff Jaureguy, McVicker Lab'
version = '1.0.0'
release = '1.0.0'

extensions = [
    'sphinx.ext.autodoc',           # Auto-generate from docstrings
    'sphinx.ext.napoleon',          # Google/NumPy docstring support
    'sphinx.ext.viewcode',          # Add source code links
    'sphinx.ext.intersphinx',       # Link to other docs
    'sphinx_autodoc_typehints',     # Use our type hints!
    'sphinx.ext.autosummary',       # Generate summary tables
]

html_theme = 'sphinx_rtd_theme'
autodoc_typehints = 'description'  # Show types in descriptions
autodoc_member_order = 'bysource'  # Keep source order
```

---

## Wave 2: Content Creation (4 Parallel Agents)

### Agent A: Basics
**Files**: `installation.rst`, `quickstart.rst`

**installation.rst** should cover:
- System requirements (bcftools, bedtools, samtools)
- Python requirements (>=3.10)
- `pip install wasp2`
- Development install
- Conda alternative

**quickstart.rst** should cover:
- 5-minute tutorial
- Basic counting example
- Basic analysis example
- Expected outputs

### Agent B: Counting Guide
**Files**: `user_guide/counting.rst`

Should cover:
- What is allele counting?
- When to use counting module
- Input requirements (BAM, VCF)
- CLI usage: `wasp2-count count-variants`
- Output format
- Common options
- Example workflows

### Agent C: Mapping Guide
**Files**: `user_guide/mapping.rst`

Should cover:
- What is WASP mapping?
- Why remapping is needed
- Input requirements
- CLI usage: `wasp2-map make-reads`, `wasp2-map filt-remapped-reads`
- Complete mapping workflow
- Common issues

### Agent D: Analysis Guide + Dev Docs
**Files**: `user_guide/analysis.rst`, `development.rst`

**analysis.rst** should cover:
- Allelic imbalance detection
- Statistical models (beta-binomial)
- CLI usage: `wasp2-analyze find-imbalance`
- Single-cell analysis
- Group comparisons
- Interpreting results

**development.rst** should cover:
- How to contribute
- Running tests
- Pre-commit hooks
- CI/CD pipeline
- Type checking with mypy

---

## Wave 3: API Documentation (3 Parallel Agents)

### Agent E: Counting API
**Files**: `api/counting.rst`

Auto-generate API docs for all counting modules:
- count_alleles.py
- count_alleles_sc.py
- filter_variant_data.py
- parse_gene_data.py
- run_counting.py
- run_counting_sc.py
- __main__.py

**Key**: Use `.. automodule::` directives to pull from type hints!

### Agent F: Mapping API
**Files**: `api/mapping.rst`

Auto-generate API docs for all mapping modules:
- filter_remap_reads.py
- intersect_variant_data.py
- make_remap_reads.py
- remap_utils.py
- run_mapping.py
- wasp_data_files.py
- __main__.py

### Agent G: Analysis API
**Files**: `api/analysis.rst`

Auto-generate API docs for all analysis modules:
- as_analysis.py (core statistical engine)
- as_analysis_sc.py
- compare_ai.py
- count_alleles.py
- count_alleles_sc.py
- filter_data.py
- run_analysis.py
- run_analysis_sc.py
- run_compare_ai.py
- __main__.py

---

## Wave 4: Integration

**Single Agent Task**

### Update index.rst
```rst
WASP2: Allele-Specific Analysis
================================

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/counting
   user_guide/mapping
   user_guide/analysis

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/counting
   api/mapping
   api/analysis

.. toctree::
   :maxdepth: 1
   :caption: Development

   development
   changelog
```

### Create .readthedocs.yaml
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
        - docs
```

### Test Build
```bash
cd docs
make clean html
# Check for warnings
# Open in browser to verify
```

---

## Critical Path

The **critical path** determines minimum time:
```
Wave 1 (foundation) → Wave 2 (longest agent) → Wave 3 (longest agent) → Wave 4 (integration)
   30 minutes            60 minutes              45 minutes               30 minutes
════════════════════════════════════════════════════════════════════════════════════════
TOTAL: 2h 45min minimum
```

All other agents run in parallel around this path.

---

## Key Success Factors

### 1. Type Hints Are Critical
Our 100% type coverage means Sphinx can auto-generate:
- Parameter types
- Return types
- Function signatures
- Type constraints

**This is why we did TH-1, TH-2, TH-3!** 🎯

### 2. Docstring Quality
If modules have good docstrings, autodoc will be beautiful.
If not, docs will still work but less descriptive.

### 3. Build Warnings
Sphinx is strict - any warnings should be fixed before committing.

---

## Agent Task Descriptions

### Wave 2 Agent A: Basics
```
Create installation.rst and quickstart.rst for WASP2 documentation.

Installation should cover:
- System deps (bcftools, bedtools, samtools)
- pip install wasp2
- conda alternative
- development install

Quickstart should have:
- 5-min counting example
- 5-min analysis example
- Expected outputs

Use .rst format with code blocks.
```

### Wave 2 Agent B: Counting
```
Create user_guide/counting.rst for WASP2 documentation.

Cover:
- Purpose of counting module
- Input requirements (BAM, VCF)
- CLI: wasp2-count count-variants
- Output format
- Common options (--samples, --region)
- Example workflow

Use .rst format with CLI examples.
```

### Wave 2 Agent C: Mapping
```
Create user_guide/mapping.rst for WASP2 documentation.

Cover:
- WASP mapping concept
- Why remapping is needed
- CLI: wasp2-map make-reads, wasp2-map filt-remapped-reads
- Complete workflow
- Common issues

Use .rst format with workflow diagrams if possible.
```

### Wave 2 Agent D: Analysis + Dev
```
Create user_guide/analysis.rst and development.rst.

Analysis should cover:
- Allelic imbalance detection
- Beta-binomial statistical model
- CLI: wasp2-analyze find-imbalance
- Single-cell analysis
- Result interpretation

Development should cover:
- Contributing guide
- Running tests
- Pre-commit hooks
- CI/CD

Use .rst format.
```

### Wave 3 Agents E/F/G: API Autodoc
```
Create api/{module}.rst with autodoc directives.

Use pattern:
.. automodule:: {module}.{file}
   :members:
   :undoc-members:
   :show-inheritance:

List all .py files in the module.
Type hints will be pulled automatically.
```

---

## Validation Checklist

After Wave 4:
- [ ] `make html` succeeds with 0 errors
- [ ] All pages render correctly
- [ ] Navigation works
- [ ] API docs show type hints
- [ ] Code examples are highlighted
- [ ] Links work (internal and external)
- [ ] Search works
- [ ] ReadTheDocs config valid

---

## Ready to Execute!

**Next steps**:
1. Execute Wave 1 (foundation) - Single agent, 30 min
2. Launch Wave 2 (content) - 4 parallel agents, 60 min
3. Launch Wave 3 (API) - 3 parallel agents, 45 min
4. Execute Wave 4 (integration) - Single agent, 30 min

**Total time**: 2h 45min - 3h 5min

Let's go! 🚀
