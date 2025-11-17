# TH-2 Implementation DAG & Parallelization Strategy

**Objective**: Reduce TH-2 implementation time from 33-42 hours to ~15-20 hours using parallel sub-agents.

---

## Dependency DAG

```
┌─────────────────────────────────────────────────────────────┐
│                    LEVEL 0 (Independent)                    │
│                    Can Run in Parallel                      │
├──────────────┬──────────────┬──────────────┬───────────────┤
│ filter_data  │count_alleles │count_alleles │  as_analysis  │
│              │              │     _sc      │  ⭐⭐⭐⭐⭐    │
│  ⭐⭐        │    ⭐⭐      │   ⭐⭐⭐      │               │
│  124 lines   │  121 lines   │  185 lines   │  476 lines    │
│  3-4 hrs     │  2-3 hrs     │  4-5 hrs     │  8-10 hrs     │
└──────┬───────┴──────┬───────┴──────┬───────┴────────┬──────┘
       │              │              │                │
       │              │              │                ▼
       │              │              │    ┌───────────────────────┐
       │              │              │    │       LEVEL 1         │
       │              │              │    │  (Depends on Level 0) │
       │              │              │    ├─────────┬─────────────┤
       │              │              └───▶│as_ana   │ compare_ai  │
       │              │                   │lysis_sc │             │
       │              │                   │  ⭐⭐⭐⭐│  ⭐⭐⭐⭐⭐   │
       │              │                   │256 lines│  512 lines  │
       │              │                   │6-8 hrs  │  10-12 hrs  │
       │              │                   └────┬────┴──────┬──────┘
       │              │                        │           │
       │              └───────────────────────▶│           │
       │                                       │           │
       │                               ┌───────▼───────┐   │
       │                               │   LEVEL 1.5   │   │
       │                               │run_analysis.py│   │
       │                               │     ⭐⭐      │   │
       │                               │   143 lines   │   │
       │                               │   2-3 hrs     │   │
       │                               └───────┬───────┘   │
       │                                       │           │
       │                                       ▼           │
       │                               ┌───────────────────▼───┐
       │                               │      LEVEL 2          │
       │                               │  run_analysis_sc.py   │
       │                               │        ⭐⭐⭐⭐        │
       │                               │      266 lines        │
       │                               │       5-6 hrs         │
       │                               └───────────┬───────────┘
       │                                           │
       │                                           ▼
       │                               ┌───────────────────────┐
       │                               │      LEVEL 3          │
       │                               │  run_compare_ai.py    │
       │                               │         ⭐⭐          │
       │                               │       73 lines        │
       │                               │       2-3 hrs         │
       │                               └───────────┬───────────┘
       │                                           │
       └───────────────────────────────────────────▼───────────┐
                                           ┌───────────────────┐
                                           │     LEVEL 4       │
                                           │   __main__.py     │
                                           │       ⭐          │
                                           │    351 lines      │
                                           │     2-3 hrs       │
                                           └───────────────────┘
```

---

## Parallelization Strategy

### Wave 1: Independent Foundation (Parallel - 4 Agents)
**Duration**: 8-10 hours (longest task in parallel set)

**Agent A1 - Quick I/O Files** (3-4 hours):
- `filter_data.py` (124 lines, ⭐⭐)
  - pysam VariantFile, BedTool operations
  - Simple I/O, no complex algorithms

**Agent A2 - Counting Module** (4-5 hours):
- `count_alleles.py` (121 lines, ⭐⭐)
- `count_alleles_sc.py` (185 lines, ⭐⭐⭐)
  - Both pileup-based counting
  - Counter and sparse array operations
  - Can be done together by one agent

**Agent A3 - Core Statistical Engine** (8-10 hours):
- `as_analysis.py` (476 lines, ⭐⭐⭐⭐⭐)
  - **MOST CRITICAL PATH** - all other files depend on this
  - scipy.optimize.OptimizeResult
  - Beta-binomial distributions
  - Likelihood ratio tests
  - Lambda functions in optimization

**Agent A4 - Wrapper** (2-3 hours):
- `run_analysis.py` (143 lines, ⭐⭐)
  - Depends on as_analysis, but lightweight wrapper
  - Can start after as_analysis signatures are known
  - **STRATEGY**: Start immediately after Agent A3 defines as_analysis function signatures

### Wave 2: Dependent Modules (Parallel - 2 Agents)
**Duration**: 10-12 hours (longest task)
**Wait for**: `as_analysis.py` completion

**Agent B1 - Single-Cell Statistical** (6-8 hours):
- `as_analysis_sc.py` (256 lines, ⭐⭐⭐⭐)
  - Imports: opt_prob, opt_phased_new, opt_unphased_dp from as_analysis
  - AnnData slicing and sparse matrices
  - Single-cell specific statistics

**Agent B2 - Group Comparison** (10-12 hours):
- `compare_ai.py` (512 lines, ⭐⭐⭐⭐⭐)
  - **LONGEST TASK IN WAVE 2**
  - Imports: opt_prob, opt_unphased_dp, opt_phased_new from as_analysis
  - Chi-square tests, FDR control
  - Complex nested returns: dict[tuple[str, str], pd.DataFrame]

### Wave 3: High-Level Runners (Sequential)
**Duration**: 5-6 hours
**Wait for**: Wave 2 completion

**Agent C1 - SC Runner** (5-6 hours):
- `run_analysis_sc.py` (266 lines, ⭐⭐⭐⭐)
  - Imports: get_imbalance_sc, adata_count_qc from as_analysis_sc
  - WaspAnalysisSC class
  - Cannot start until as_analysis_sc is complete

### Wave 4: Integration Layer (Sequential)
**Duration**: 2-3 hours
**Wait for**: Wave 3 completion

**Agent D1 - Comparison Runner** (2-3 hours):
- `run_compare_ai.py` (73 lines, ⭐⭐)
  - Imports: adata_count_qc, process_adata_inputs, get_compared_imbalance
  - Small wrapper - just wiring components together

### Wave 5: CLI (Sequential)
**Duration**: 2-3 hours
**Wait for**: All previous waves

**Agent E1 - CLI** (2-3 hours):
- `__main__.py` (351 lines, ⭐)
  - Already 95% typed (Typer does most work)
  - Just needs wrapper returns typed

---

## Time Analysis

### Sequential Implementation (Original Plan)
```
Tier 1: 2-3 hrs
Tier 2: 3-4 hrs
Tier 3: 4-5 hrs + 5-6 hrs = 9-11 hrs
Tier 4: 6-8 hrs
Tier 5: 8-10 hrs
Tier 6: 10-12 hrs
───────────────
TOTAL: 42-54 hours
```

### Parallel Implementation (Optimized)
```
Wave 1 (parallel, 4 agents): max(3-4, 4-5, 8-10, 2-3) = 8-10 hrs
Wave 2 (parallel, 2 agents): max(6-8, 10-12) = 10-12 hrs
Wave 3 (sequential):          5-6 hrs
Wave 4 (sequential):          2-3 hrs
Wave 5 (sequential):          2-3 hrs
────────────────────────────────────────────────────────
TOTAL: 27-34 hours
```

**Time Savings**: 8-15 hours (21-37% reduction)

---

## Critical Path Analysis

**Critical Path** (determines minimum time):
```
as_analysis.py → compare_ai.py → run_compare_ai.py → __main__.py
   8-10 hrs       10-12 hrs         2-3 hrs          2-3 hrs
═══════════════════════════════════════════════════════════
TOTAL CRITICAL PATH: 22-28 hours
```

All other files can be worked in parallel around this path.

---

## Sub-Agent Implementation Plan

### Phase 1: Launch Wave 1 (4 parallel agents)

**Agent A1**: "Quick I/O Files"
```
Task: Add type hints to filter_data.py
- pysam.VariantFile operations
- Path handling
- Simple I/O patterns
Duration: 3-4 hours
```

**Agent A2**: "Counting Module"
```
Task: Add type hints to count_alleles.py and count_alleles_sc.py
- pileup operations
- Counter typing
- SparseArray from pandas
Duration: 4-5 hours
```

**Agent A3**: "Core Statistical Engine" ⚠️ CRITICAL PATH
```
Task: Add type hints to as_analysis.py
- scipy.optimize.OptimizeResult
- Beta-binomial distributions
- Likelihood ratio tests
- Lambda function typing
- NDArray[np.float64]
Duration: 8-10 hours
Priority: HIGHEST - blocks Wave 2
```

**Agent A4**: "Wrapper"
```
Task: Add type hints to run_analysis.py
- Wait for as_analysis.py function signatures (30 min into A3)
- Lightweight wrapper functions
- Polars DataFrame operations
Duration: 2-3 hours
Dependencies: Partial on A3 (just needs signatures)
```

### Phase 2: Launch Wave 2 (2 parallel agents)
**Trigger**: As soon as Agent A3 (as_analysis.py) completes

**Agent B1**: "Single-Cell Statistical"
```
Task: Add type hints to as_analysis_sc.py
- Uses: opt_prob, opt_phased_new, opt_unphased_dp from as_analysis
- AnnData slicing
- Sparse matrix operations
Duration: 6-8 hours
```

**Agent B2**: "Group Comparison" ⚠️ CRITICAL PATH
```
Task: Add type hints to compare_ai.py
- Uses: opt_prob, opt_unphased_dp, opt_phased_new from as_analysis
- Chi-square tests
- FDR control
- Complex nested dict returns
Duration: 10-12 hours
Priority: HIGHEST - blocks Wave 3
```

### Phase 3: Launch Wave 3 (1 agent)
**Trigger**: Agent B1 (as_analysis_sc.py) completes

**Agent C1**: "SC Runner"
```
Task: Add type hints to run_analysis_sc.py
- Uses: get_imbalance_sc, adata_count_qc from as_analysis_sc
- WaspAnalysisSC class
Duration: 5-6 hours
```

### Phase 4: Launch Wave 4 (1 agent)
**Trigger**: Agent C1 (run_analysis_sc.py) completes

**Agent D1**: "Comparison Runner"
```
Task: Add type hints to run_compare_ai.py
- Uses: adata_count_qc, process_adata_inputs, get_compared_imbalance
- Small wrapper file
Duration: 2-3 hours
```

### Phase 5: Launch Wave 5 (1 agent)
**Trigger**: All Waves 1-4 complete

**Agent E1**: "CLI"
```
Task: Add type hints to __main__.py
- Already 95% typed (Typer handles most)
- Just wrapper return types
Duration: 2-3 hours
```

---

## Risk Mitigation

### Risk 1: Agent A3 (as_analysis.py) Delays
**Impact**: Blocks entire Wave 2 (16-20 hours of work)
**Mitigation**:
- Start Agent A3 FIRST in Wave 1
- Monitor progress closely
- Have backup plan to split file if needed
- Prioritize getting function signatures defined early

### Risk 2: Type Inconsistencies Between Waves
**Impact**: Wave 2 agents use wrong types from Wave 1
**Mitigation**:
- After Wave 1 completes, run `mypy src/analysis/` to validate
- Share type definitions file between agents
- Use explicit imports in task descriptions

### Risk 3: Agent B2 (compare_ai.py) Complexity
**Impact**: Longest single file (10-12 hours)
**Mitigation**:
- Consider splitting into two sub-tasks if agent struggles
- Provide detailed examples of complex return types
- Reference as_analysis.py types explicitly

---

## Validation Strategy

After each wave:
```bash
# Type check completed files
mypy src/analysis/filter_data.py src/analysis/count_alleles.py ...

# Run tests to ensure no runtime breakage
python -m pytest tests/ -v

# Commit working set
git add src/analysis/*.py
git commit -m "TH-2: Complete Wave N - [file list]"
```

---

## Expected Outcome

**Before TH-2**: 14/24 files typed (58% complete)
**After TH-2**: 24/24 files typed (100% complete) ✅

**Time Comparison**:
- Sequential: 42-54 hours
- Parallel (5 waves): 27-34 hours
- **Savings: 8-20 hours (21-37% reduction)**

**Quality**: All files pass mypy validation, all tests passing

---

## Recommended Execution

**Option 1: Aggressive Parallel (4 agents Wave 1)**
- Fastest overall time (27-34 hours)
- Requires monitoring 4 concurrent agents in Wave 1
- Best if you can manage multiple agent outputs

**Option 2: Conservative Parallel (2 agents Wave 1)**
- Launch only A2 (counting) + A3 (core engine) in parallel
- Do A1 and A4 sequentially after
- Extends Wave 1 to 11-15 hours, but easier to manage
- Total time: 31-38 hours

**Option 3: Critical Path Focus (1 agent + helpers)**
- Start with A3 (as_analysis.py) immediately
- Launch other Wave 1 agents only after A3 is 50% complete
- Safest approach, ensures critical path completes
- Total time: 35-42 hours

**Recommendation**: **Option 1** - You've successfully managed complex multi-file implementations before. The 8-20 hour time savings is worth the coordination overhead.

---

## Ready to Execute?

**Next Steps**:
1. Review DAG and strategy
2. Choose execution option (1, 2, or 3)
3. Launch Wave 1 agents with detailed task descriptions
4. Monitor critical path (Agent A3)
5. Validate after each wave
6. Celebrate 100% type coverage! 🎉
