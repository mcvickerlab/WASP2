# TH-2 Agent Launch Guide

Quick reference for launching parallel sub-agents for TH-2 implementation.

---

## Wave 1: Launch All 4 Agents Simultaneously

### Agent A1 - Quick I/O Files (3-4 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/filter_data.py (124 lines).

Focus areas:
- pysam.VariantFile → VariantFile type
- pysam.AlignmentFile → AlignmentFile type
- pybedtools.BedTool → BedTool type
- Path operations → Union[str, Path]
- Function return types → str for file paths

Key patterns from TH-1 and TH-3:
- Use Optional[str] for optional parameters
- AlignmentFile for BAM operations
- VariantFile for VCF operations

After completion:
- Run: mypy src/analysis/filter_data.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

### Agent A2 - Counting Module (4-5 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/count_alleles.py (121 lines)
and src/analysis/count_alleles_sc.py (185 lines).

Focus areas:
- pysam pileup operations → pileup_col typing
- collections.Counter → Counter[str]
- pandas.arrays.SparseArray → SparseArray type
- Tuple returns → Tuple[List[str], List[str]]
- AlignmentFile parameters

Key patterns:
count_alleles.py:
- pileup_pos() returns Tuple[List[str], List[str]]
- AlignmentFile from pysam.libcalignmentfile

count_alleles_sc.py:
- parse_barcode() → Optional[str] (can return None)
- SparseArray from pandas.arrays
- np.ndarray for numpy arrays

After completion:
- Run: mypy src/analysis/count_alleles.py src/analysis/count_alleles_sc.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

### Agent A3 - Core Statistical Engine ⚠️ CRITICAL (8-10 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/as_analysis.py (476 lines).

⚠️ CRITICAL PATH ⚠️ - All Wave 2 files depend on this completion.

Focus areas:
1. scipy.optimize returns → scipy.optimize.OptimizeResult
2. Beta-binomial distributions → scipy.stats.betabinom type
3. Likelihood ratio tests → float returns
4. Lambda functions → Callable[[float], float] or Callable[..., float]
5. Numpy arrays → NDArray[np.float64]

Key functions to type:
- opt_prob() - optimization function
- opt_phased_new() - phased variant optimization
- opt_unphased_dp() - unphased variant optimization
- get_imbalance() - main analysis function

Critical patterns:
```python
from scipy.optimize import minimize, minimize_scalar, OptimizeResult
from scipy.stats import betabinom, chi2
from numpy.typing import NDArray
import numpy as np

def opt_prob(...) -> OptimizeResult:
    result = minimize(...)
    return result

def likelihood_function(params: NDArray[np.float64]) -> float:
    ...
    return -log_likelihood

# Lambda typing
optimizer: Callable[[NDArray[np.float64]], float] = lambda x: -log_likelihood(x)
```

After completion:
- Run: mypy src/analysis/as_analysis.py
- Verify: 0 errors (may have warnings about scipy stubs - ignore those)
- Document function signatures for Wave 2 agents
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

### Agent A4 - Wrapper (2-3 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/run_analysis.py (143 lines).

Dependencies:
- Wait 30 minutes for Agent A3 to define as_analysis.py function signatures
- Then proceed with run_analysis.py typing

Focus areas:
- Imports from as_analysis → get_imbalance function
- polars.DataFrame operations → pl.DataFrame type
- Path and file I/O → Union[str, Path]
- Wrapper function returns → polars.DataFrame or None

Key patterns:
```python
import polars as pl
from as_analysis import get_imbalance

def run_ai_analysis(...) -> Optional[pl.DataFrame]:
    # Process data
    result = get_imbalance(...)  # Use types from as_analysis.py
    return result
```

After completion:
- Run: mypy src/analysis/run_analysis.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

## Wave 2: Launch After Agent A3 Completes

### Agent B1 - Single-Cell Statistical (6-8 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/as_analysis_sc.py (256 lines).

Prerequisites:
- Agent A3 (as_analysis.py) MUST be complete
- Import signatures from as_analysis: opt_prob, opt_phased_new, opt_unphased_dp

Focus areas:
- AnnData operations → anndata.AnnData type
- Sparse matrices → scipy.sparse matrices
- Single-cell specific statistics
- Uses optimization functions from as_analysis

Key patterns:
```python
import anndata
from scipy import sparse
from as_analysis import opt_prob, opt_phased_new, opt_unphased_dp

def get_imbalance_sc(adata: anndata.AnnData, ...) -> anndata.AnnData:
    # Use opt_prob, opt_phased_new from as_analysis
    result = opt_prob(...)  # Already typed in Wave 1
    return adata
```

After completion:
- Run: mypy src/analysis/as_analysis_sc.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

### Agent B2 - Group Comparison ⚠️ CRITICAL (10-12 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/compare_ai.py (512 lines).

Prerequisites:
- Agent A3 (as_analysis.py) MUST be complete
- Import signatures from as_analysis: opt_prob, opt_unphased_dp, opt_phased_new

⚠️ LONGEST SINGLE FILE ⚠️ - Most complex return types in the module.

Focus areas:
1. Chi-square tests → scipy.stats.chi2 typing
2. False discovery rate → scipy.stats.false_discovery_control
3. Complex nested returns → dict[tuple[str, str], pd.DataFrame]
4. Group comparison statistics
5. Uses all three opt_* functions from as_analysis

Key patterns:
```python
from scipy.stats import chi2, false_discovery_control
from as_analysis import opt_prob, opt_unphased_dp, opt_phased_new
import pandas as pd

def get_compared_imbalance(...) -> dict[tuple[str, str], pd.DataFrame]:
    """
    Returns dict mapping (group1, group2) tuples to comparison DataFrames
    """
    results: dict[tuple[str, str], pd.DataFrame] = {}
    # ... analysis
    return results
```

After completion:
- Run: mypy src/analysis/compare_ai.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

## Wave 3: Launch After Agent B1 Completes

### Agent C1 - SC Runner (5-6 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/run_analysis_sc.py (266 lines).

Prerequisites:
- Agent B1 (as_analysis_sc.py) MUST be complete
- Import signatures: get_imbalance_sc, adata_count_qc from as_analysis_sc

Focus areas:
- WaspAnalysisSC class typing
- AnnData pipeline operations
- process_adata_inputs() function
- Class method return types

Key patterns:
```python
from as_analysis_sc import get_imbalance_sc, adata_count_qc
import anndata

class WaspAnalysisSC:
    def __init__(self, ...) -> None:
        ...

    def run_analysis(self) -> anndata.AnnData:
        result = get_imbalance_sc(...)  # Already typed in Wave 2
        return result

def process_adata_inputs(...) -> anndata.AnnData:
    ...
```

After completion:
- Run: mypy src/analysis/run_analysis_sc.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

## Wave 4: Launch After Agent C1 Completes

### Agent D1 - Comparison Runner (2-3 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/run_compare_ai.py (73 lines).

Prerequisites:
- Agent B1 (as_analysis_sc.py) complete
- Agent B2 (compare_ai.py) complete
- Agent C1 (run_analysis_sc.py) complete

Focus areas:
- Small wrapper file - just wiring components
- Imports: adata_count_qc, process_adata_inputs, get_compared_imbalance
- Simple function returns

Key patterns:
```python
from as_analysis_sc import adata_count_qc
from run_analysis_sc import process_adata_inputs
from compare_ai import get_compared_imbalance

def run_ai_comparison(...) -> dict[tuple[str, str], pd.DataFrame]:
    # Wire together components
    result = get_compared_imbalance(...)
    return result
```

After completion:
- Run: mypy src/analysis/run_compare_ai.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

## Wave 5: Launch After All Waves 1-4 Complete

### Agent E1 - CLI (2-3 hours)

**Task Description**:
```
Add comprehensive type hints to src/analysis/__main__.py (351 lines).

Prerequisites:
- ALL previous waves complete
- Imports from: run_analysis, run_analysis_sc, run_compare_ai

Focus areas:
- Typer CLI commands (already 95% typed)
- Wrapper function return types
- Command decorators (Typer handles most typing)

Key patterns:
```python
import typer
from typing import Annotated, Optional
from run_analysis import run_ai_analysis
from run_analysis_sc import run_ai_analysis_sc
from run_compare_ai import run_ai_comparison

app = typer.Typer()

@app.command()
def analyze(...) -> None:
    """CLI command - returns None"""
    run_ai_analysis(...)
```

After completion:
- Run: mypy src/analysis/__main__.py
- Verify: 0 errors
- Run: python -m pytest tests/ -v
- Verify: All tests pass
```

---

## Final Validation

After all waves complete:

```bash
# Type check entire analysis module
mypy src/analysis/

# Should show: Success: no issues found in 10 source files

# Run full test suite
python -m pytest tests/ -v

# Should show: All tests passing

# Check type coverage across entire project
mypy src/

# Should show: Success: no issues found in 24 source files

# Commit all changes
git add src/analysis/*.py
git commit -m "TH-2: Add comprehensive type hints to analysis module (100% coverage)"
git push -u origin claude/merge-test-data-bundle-01SeJo12Zuj6GVrFoQxw9HiG
```

---

## Success Criteria

✅ All 10 analysis files fully typed
✅ mypy src/analysis/ returns 0 errors
✅ All regression tests passing
✅ 100% type coverage: 24/24 files complete
✅ Estimated time: 27-34 hours (vs 42-54 sequential)
✅ Time savings: 8-20 hours (21-37% reduction)

---

## Notes

- **Critical Path**: A3 → B2 determines minimum completion time (22-28 hours)
- **Monitor**: Agent A3 and B2 closely - they're the longest tasks
- **Validate**: Run mypy after each wave to catch type mismatches early
- **Commit**: After each wave completion to save progress
- **Celebrate**: 100% type coverage is a major milestone! 🎉
