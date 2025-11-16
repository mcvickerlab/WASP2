# Regression Test Suite

**Purpose:** Validate that code changes don't break functionality or degrade performance.

## Quick Start

```bash
# Run all regression tests
pytest tests/regression/ -v

# Run specific test class
pytest tests/regression/test_pipeline_regression.py::TestCountingRegression -v

# Run with performance tests (slow)
pytest tests/regression/ -v -m slow
```

## What Gets Tested

### ✅ Output Correctness
- **MD5 checksums** - Outputs must match baseline exactly
- **File structure** - Column names, data types, row counts
- **Statistical validity** - Values in correct ranges (p-values [0,1], etc.)

### ⚡ Performance
- **Memory usage** - Must not exceed baseline × 1.20 (20% tolerance)
- **Execution time** - Must not exceed baseline × 1.30 (30% tolerance)
- **WASP filter rate** - Must keep >95% of reads

### 📊 Baselines Used

From `baselines/` directory (committed):
```
Counting:  9.26s, 639 MB, MD5: 127a81810a43db3cc6924a26f591cc7a
Analysis:  2.97s, 340 MB, MD5: 394e1a7dbf14220079c3142c5b15bad8
Mapping:   8s,    488 MB, 125,387 reads kept (99%)
```

## Usage Workflow

### Before Refactoring
```bash
# Ensure all tests pass
pytest tests/regression/ -v

# If any fail, investigate before starting
```

### During Refactoring
```bash
# Run tests frequently (after each logical change)
pytest tests/regression/ -v

# Run fast tests only (skip full pipeline)
pytest tests/regression/ -v -m "not slow"
```

### After Refactoring
```bash
# Run full test suite including slow E2E tests
pytest tests/regression/ -v -m slow

# If MD5 changed but output is correct, update baseline:
# 1. Manually verify new output is correct
# 2. Update MD5 in test_pipeline_regression.py:BASELINE_EXPECTATIONS
# 3. Commit new baseline files
```

## Test Categories

| Test Class | Speed | What It Tests |
|------------|-------|---------------|
| `TestCountingRegression` | Fast (1s) | Counting output, memory, performance |
| `TestAnalysisRegression` | Fast (1s) | Analysis output, memory, performance |
| `TestMappingRegression` | Fast (1s) | WASP filtering, read counts |
| `TestFullPipelineIntegration` | Slow (20s) | End-to-end reproducibility |

## Continuous Integration

Add to `.github/workflows/regression.yml`:

```yaml
name: Regression Tests
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.11'
      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest pandas
      - name: Run regression tests
        run: pytest tests/regression/ -v
```

## Updating Baselines

When you **intentionally** change outputs:

1. **Verify change is correct**
   ```bash
   # Compare old vs new output
   diff baselines/counting/counts.tsv new_output/counts.tsv
   ```

2. **Update baseline files**
   ```bash
   # Run pipeline to regenerate baselines
   ./scripts/run_full_pipeline_baseline.sh
   ```

3. **Update expected MD5s**
   ```bash
   # Calculate new checksums
   md5sum baselines/counting/counts.tsv
   md5sum baselines/analysis/ai_results.tsv

   # Update BASELINE_EXPECTATIONS in test_pipeline_regression.py
   ```

4. **Commit changes**
   ```bash
   git add baselines/ tests/regression/test_pipeline_regression.py
   git commit -m "Update baselines after [description of change]"
   ```

## Troubleshooting

### Test fails with MD5 mismatch
**Cause:** Output has changed
**Fix:** Compare outputs to verify correctness, then update baseline

### Test fails with memory regression
**Cause:** Code now uses more memory
**Fix:** Investigate memory leak or optimize, OR increase tolerance if justified

### Test fails with performance regression
**Cause:** Code is slower
**Fix:** Profile and optimize hot paths, OR increase tolerance if complexity trade-off

### Test skipped
**Cause:** Baseline files not found
**Fix:** Run `./scripts/run_full_pipeline_baseline.sh` to generate baselines

## Philosophy

> **"Tests are a safety net, not a straightjacket"**

- ✅ Tests should **enable** refactoring, not prevent it
- ✅ Tolerances exist to avoid flaky tests (±20-30%)
- ✅ Update baselines when outputs **intentionally** change
- ❌ Don't disable tests just because they fail
- ❌ Don't increase tolerances to paper over problems

## See Also

- `baselines/pipeline_metadata.txt` - Detailed benchmark data
- `docs/modules/COUNTING_MODULE.md` - Module documentation
- `ENGINEERING_PLAN.md` - Overall refactoring strategy
