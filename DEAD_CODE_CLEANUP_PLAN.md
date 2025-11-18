# WASP2 Dead Code Cleanup Plan

**Date:** 2025-11-18
**Total Dead Code Found:** ~900+ lines
**Risk Assessment:** LOW (all changes are safe)

---

## 📊 Summary Statistics

### Files Affected
- **3 entire legacy files** to DELETE (438 lines)
- **7 large commented blocks** to REMOVE (300+ lines)
- **4 unused functions** to DELETE (180+ lines)
- **2 unused imports** to REMOVE

### By Module
| Module | Dead Code Lines | Files Affected |
|--------|----------------|----------------|
| Mapping | ~200 | 3 files |
| Counting | ~0 | 0 files (clean!) |
| Analysis | ~700 | 7 files |

---

## 🎯 PHASE 1: IMMEDIATE CLEANUP (SAFE)

### Priority 1: Delete Entire Legacy Files

**ANALYSIS MODULE - 3 Files to DELETE:**

1. **`src/analysis/filter_data.py`** (126 lines)
   - Status: LEGACY from Python 3.8
   - Reason: Superseded by counting module
   - Risk: LOW - No imports found in codebase
   - Action: `git rm src/analysis/filter_data.py`

2. **`src/analysis/count_alleles.py`** (124 lines)
   - Status: LEGACY from Python 3.8
   - Reason: Superseded by counting module
   - Risk: LOW - No imports found in codebase
   - Action: `git rm src/analysis/count_alleles.py`

3. **`src/analysis/count_alleles_sc.py`** (188 lines)
   - Status: LEGACY from Python 3.8
   - Reason: Superseded by counting module
   - Risk: LOW - No imports found in codebase
   - Action: `git rm src/analysis/count_alleles_sc.py`

**Validation:**
```bash
# Verify no imports exist
grep -r "from analysis.filter_data" src/
grep -r "from analysis.count_alleles" src/
grep -r "import filter_data" src/
grep -r "import count_alleles" src/
# Should return: no matches
```

---

### Priority 2: Remove Unused Imports

**File: `src/mapping/filter_remap_reads.py`**
- Line 3: `import timeit`
- Reason: Never used in file
- Action: Delete line 3

**File: `src/analysis/count_alleles_sc.py`** *(will be deleted in Priority 1)*
- Line 16: `from pysam import VariantFile`
- Reason: Only in comment
- Action: N/A (file deleted)

---

### Priority 3: Remove Unused Functions

**File: `src/analysis/as_analysis.py`**

1. **Function: `opt_phased()`** (Lines 81-104, 24 lines)
   - Reason: Replaced by `opt_phased_new()` at line 110
   - Risk: LOW - New version is used
   - Action: Delete lines 81-104

2. **Function: `opt_unphased()`** (Lines 137-160, 24 lines)
   - Reason: Replaced by `opt_unphased_dp()` at line 164
   - Risk: LOW - New version is used
   - Action: Delete lines 137-160

3. **Function: `get_imbalance_sc()`** (Lines 482-602, 121 lines)
   - Reason: Marked as "LEGACY, NOT REALLY USED"
   - Replaced by: `as_analysis_sc.py` version
   - Risk: LOW - Modern version exists
   - Action: Delete lines 482-602

**File: `src/analysis/compare_ai.py`**

4. **Function: `get_imbalance_func()`** (Lines 23-53, 31 lines)
   - Reason: Only called by unused V0 function
   - Risk: LOW - Not in active code path
   - Action: Delete lines 23-53 (will remove with V0 functions)

---

### Priority 4: Remove Large Commented Code Blocks

**File: `src/mapping/make_remap_reads.py`**

1. **Lines 422-525 (103 lines)** - Old `swap_chrom_alleles()` implementation
   - Reason: Active implementation exists at line 135
   - Action: Delete commented block

**File: `src/mapping/intersect_variant_data.py`**

2. **Lines 108-141 (34 lines)** - Old `process_bam()` implementation
   - Reason: Active version at line 71
   - Action: Delete commented block

3. **Lines 161-184 (24 lines)** - Old `filter_intersect_data()` wrapper
   - Reason: Logic moved to run_mapping.py
   - Action: Delete commented block

4. **Lines 188-244 (57 lines)** - Old `make_intersect_df()` implementation
   - Reason: Optimized version at line 247
   - Action: Delete commented block

**File: `src/mapping/remap_utils.py`**

5. **Lines 109-122 (14 lines)** - Old `get_read_het_data()` signature
   - Reason: Updated version at line 88
   - Action: Delete commented block

**File: `src/mapping/make_remap_reads.py`**

6. **Line 9** - Single line: `# from collections import defaultdict`
   - Action: Delete

---

## 🔍 PHASE 2: REVIEW BEFORE REMOVAL (DECISION NEEDED)

### V0 Comparison Functions (230 lines)

**File: `src/analysis/compare_ai.py`**

Two V0 functions marked "COULD BE USEFUL AS AN OPTION POSSIBLY":

1. **`get_compared_imbalance_diff_snps()`** (Lines 368-472, 105 lines)
   - What: V0 version without shared SNPs requirement
   - Used: Never called
   - Decision Needed: Keep as alternative approach or remove?
   - Recommendation: **REMOVE** - If needed later, recover from git history

2. **`compare_imbalance_between_groups_diff_snps()`** (Lines 475-595, 121 lines)
   - What: Helper for above function
   - Used: Only by function #1
   - Decision: Remove with parent function

**Action Plan:**
- [ ] Confirm with team: Is "diff SNPs" approach needed?
- [ ] If NO: Remove both functions (230 lines)
- [ ] If YES: Refactor into proper optional parameter in main function

---

## ✅ PHASE 3: TESTING & VALIDATION PLAN

### Pre-Cleanup Baseline

**Step 1: Run Full Test Suite**
```bash
# Activate environment
conda activate WASP2

# Run all tests
pytest tests/ -v

# Run regression tests
pytest tests/regression/ -v

# Run type checking
mypy src/

# Record results
```

**Step 2: Run Full Pipeline Baseline**
```bash
# Run complete pipeline with test data
bash scripts/run_full_pipeline_baseline.sh

# Save outputs as baseline
cp baselines/counting/counts_head20.txt /tmp/baseline_counts_before.txt
cp baselines/analysis/ai_results_head20.txt /tmp/baseline_ai_before.txt
```

---

### During Cleanup: Incremental Validation

**For Each Phase:**

1. **Make changes**
2. **Run type check immediately:**
   ```bash
   mypy src/
   # Should: 0 errors
   ```

3. **Run quick sanity tests:**
   ```bash
   # Test imports work
   python -c "from counting import count_alleles"
   python -c "from mapping import remap_utils"
   python -c "from analysis import as_analysis"
   ```

4. **Run affected module tests:**
   ```bash
   # If mapping changed:
   pytest tests/regression/test_pipeline_regression.py::TestMapping -v

   # If analysis changed:
   pytest tests/regression/test_pipeline_regression.py::TestAnalysis -v
   ```

---

### Post-Cleanup Full Validation

**Step 1: Complete Test Suite**
```bash
# All tests should still pass
pytest tests/ -v

# Regression tests
pytest tests/regression/ -v

# Type checking
mypy src/

# All should: PASS with same results as pre-cleanup
```

**Step 2: Full Pipeline Verification**
```bash
# Re-run full pipeline
bash scripts/run_full_pipeline_baseline.sh

# Compare outputs byte-for-byte
diff /tmp/baseline_counts_before.txt baselines/counting/counts_head20.txt
diff /tmp/baseline_ai_before.txt baselines/analysis/ai_results_head20.txt

# Should: No differences (identical outputs)
```

**Step 3: Manual CLI Testing**
```bash
# Test each CLI entry point
python -m counting count-variants --help
python -m mapping make-reads --help
python -m analysis find-imbalance --help

# Should: All display help correctly
```

**Step 4: Import Testing**
```bash
# Test all public APIs still work
python -c "
from counting.count_alleles import count_snp_alleles
from mapping.remap_utils import make_phased_seqs
from analysis.as_analysis import get_imbalance
print('All imports successful!')
"
```

---

## 🚀 EXECUTION PLAN

### Timeline: 1-2 hours

**Phase 1A: Delete Legacy Files (10 min)**
```bash
# Create cleanup branch
git checkout -b cleanup/remove-dead-code

# Delete 3 legacy files
git rm src/analysis/filter_data.py
git rm src/analysis/count_alleles.py
git rm src/analysis/count_alleles_sc.py

# Test
mypy src/analysis/
pytest tests/regression/ -v

# Commit
git commit -m "Cleanup: Remove 3 legacy analysis files (438 lines)"
```

**Phase 1B: Remove Unused Imports (5 min)**
```bash
# Edit files manually or with sed
sed -i '3d' src/mapping/filter_remap_reads.py  # Delete line 3

# Test
mypy src/mapping/
pytest tests/regression/test_pipeline_regression.py::TestMapping -v

# Commit
git commit -am "Cleanup: Remove unused import in filter_remap_reads"
```

**Phase 1C: Remove Unused Functions (15 min)**
```bash
# Edit src/analysis/as_analysis.py
# Delete lines 81-104, 137-160, 482-602

# Edit src/analysis/compare_ai.py
# Delete lines 23-53 (if removing V0 functions)

# Test
mypy src/analysis/
pytest tests/regression/test_pipeline_regression.py::TestAnalysis -v

# Commit
git commit -am "Cleanup: Remove unused optimization functions (180 lines)"
```

**Phase 1D: Remove Commented Blocks (15 min)**
```bash
# Edit each file, remove commented blocks:
# - src/mapping/make_remap_reads.py: lines 422-525, line 9
# - src/mapping/intersect_variant_data.py: lines 108-141, 161-184, 188-244
# - src/mapping/remap_utils.py: lines 109-122

# Test
mypy src/mapping/
pytest tests/regression/test_pipeline_regression.py::TestMapping -v

# Commit
git commit -am "Cleanup: Remove commented dead code blocks (300 lines)"
```

**Phase 2: Full Validation (20 min)**
```bash
# Run complete test suite
pytest tests/ -v

# Run full pipeline
bash scripts/run_full_pipeline_baseline.sh

# Compare against baseline
bash scripts/validate_against_baseline.sh

# If all pass: Ready to merge!
```

**Phase 3: Merge & Push (5 min)**
```bash
# Push cleanup branch
git push -u origin cleanup/remove-dead-code

# Merge to main development branch
git checkout claude/explore-codebase-01XDRjqauxDuSFC3nPBdG4P3
git merge cleanup/remove-dead-code

# Push
git push
```

---

## 📋 SANITY CHECKS

### Before Each Commit:
- [ ] mypy src/ → 0 errors
- [ ] Affected tests pass
- [ ] Code still runs

### Before Final Merge:
- [ ] All tests pass (pytest tests/ -v)
- [ ] Type checking passes (mypy src/)
- [ ] Full pipeline runs successfully
- [ ] Outputs match baseline exactly
- [ ] All CLI commands work (--help)
- [ ] All public imports work

### Rollback Plan:
If anything breaks:
```bash
# Revert to previous commit
git reset --hard HEAD~1

# Or revert specific file
git checkout HEAD~1 -- src/path/to/file.py
```

---

## 📈 EXPECTED RESULTS

**Lines of Code Removed:**
- ~900+ lines of dead code
- ~19 TODO comments remain (flagged for future work)
- 0 functional changes

**Files Changed:**
- 3 files deleted
- 7 files modified (comments/dead code removed)
- 0 files added

**Performance Impact:**
- Slightly faster imports (less code to parse)
- Cleaner diffs for future changes
- Easier code review

**Maintenance Impact:**
- Less confusion for new contributors
- Clearer which functions are active
- Easier to find relevant code

---

## 🎯 SUCCESS CRITERIA

✅ All tests pass (same results as before cleanup)
✅ Pipeline outputs identical to baseline
✅ mypy reports 0 errors
✅ All CLI commands functional
✅ ~900 lines removed
✅ No functional changes
✅ Clean git history with clear commit messages

---

## 📝 NOTES

- Keep TODO comments → future feature flags
- V0 functions → DECIDE: remove or keep as optional
- All changes are LOW RISK
- Full git history available if code needed later
- Consider adding test coverage for flagged areas

**Ready to execute?** Start with Phase 1A!
