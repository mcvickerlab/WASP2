# Quick Start: Statistical Analysis

**TL;DR**: Run simulation, analyze results, get publication-ready stats.

## One-Line Usage

```bash
python analyze_simulation_results.py simulation_results.csv -o report.md
```

## Complete Workflow

### Step 1: Run Simulation

```bash
# Run moderate tier (270 tests, ~30 min)
python simulate_indel_ase_v2.py --tier moderate
```

This creates: `simulation_results.csv`

### Step 2: Analyze Results

```bash
# Generate statistical report
python analyze_simulation_results.py simulation_results.csv \
    --output statistical_report.md \
    --plots ./figures/
```

This creates:
- `statistical_report.md` - Full statistical analysis
- `./figures/*.png` - Diagnostic plots

### Step 3: Use in Manuscript

Open `statistical_report.md` and copy:

1. **For Results section**: Copy "Manuscript-Ready Text > Results Section"
2. **For Methods section**: Copy "Manuscript-Ready Text > Methods Section"
3. **For Supplement**: Copy entire summary table

## What You Get

### Summary Statistics Table

| Category | N | Mean Error (%) | 95% CI | Pass Rate |
|----------|---|----------------|--------|-----------|
| Overall | 270 | 2.7 | [2.5, 3.0] | 100% |
| SNP | 90 | 2.1 | [1.9, 2.4] | 100% |
| INS | 90 | 2.5 | [2.2, 2.7] | 100% |
| DEL | 90 | 2.5 | [2.3, 2.8] | 100% |

### Hypothesis Tests

**Test 1: Bias Detection**
- H0: Mean error = 0
- Result: Small bias (2.7%) detected but acceptable

**Test 2: Threshold Validation**
- H0: Mean error ≥ 10% (FAIL)
- Result: PASS (p < 0.001)

### Consistency Metrics

- **CV**: 0.48 (low variability)
- **ICC**: 0.82 (excellent reproducibility)

### Sample Size Justification

- N=10 replicates sufficient for deterministic algorithm
- Can detect effect sizes ≥0.72 with 80% power
- Following van de Geijn et al. (2015) precedent

## Customization Options

### Change Significance Level

```bash
python analyze_simulation_results.py results.csv --alpha 0.01
```

### More Bootstrap Iterations

```bash
python analyze_simulation_results.py results.csv --bootstrap 20000
```

### Just Print to Screen

```bash
python analyze_simulation_results.py results.csv
```

## For Reviewers

If reviewers ask:

**"Why only 10 replicates?"**
→ See "Sample Size Justification" section in report

**"Is the bias significant?"**
→ Yes (p<0.001) but biologically negligible (2.7% vs 10% threshold)

**"What about multiple testing?"**
→ Primary hypothesis tested once; stratified = descriptive only

## Files Overview

- `analyze_simulation_results.py` - Main script (run this)
- `README_STATISTICAL_ANALYSIS.md` - Full documentation
- `EXAMPLE_STATISTICAL_OUTPUT.md` - Example output
- `STATISTICAL_ANALYSIS_SUMMARY.md` - Feature summary

## Help

```bash
python analyze_simulation_results.py --help
```

For detailed docs: See `README_STATISTICAL_ANALYSIS.md`

---

**Status**: Production-ready
**Branch**: `feat/simulation-statistics`
**Tested**: ✅ Works with real and synthetic data
