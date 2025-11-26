# Statistical Analysis of WASP2 Simulation Results

This document explains how to perform publication-ready statistical analysis on WASP2 simulation validation results.

## Overview

The `analyze_simulation_results.py` script provides rigorous statistical validation including:

1. **Bootstrap Confidence Intervals** (95% CI) for all error metrics
2. **Hypothesis Tests** to prove algorithmic correctness
3. **Sample Size Justification** via power analysis
4. **Variance Metrics** to demonstrate consistency
5. **Publication-Quality Output** ready for manuscript inclusion

## Quick Start

### Basic Usage

```bash
# Run simulation first
python simulate_indel_ase_v2.py --tier moderate

# Analyze results
python analyze_simulation_results.py simulation_results.csv
```

This will print a complete statistical report to stdout.

### Save Report to File

```bash
python analyze_simulation_results.py simulation_results.csv --output analysis_report.md
```

### Generate Diagnostic Plots

```bash
python analyze_simulation_results.py simulation_results.csv --plots ./figures/
```

This creates:
- `error_distribution.png` - Histogram showing error distribution
- `qq_plot.png` - Q-Q plot for normality assessment
- `error_by_variant_type.png` - Boxplot comparing SNP/INS/DEL
- `error_by_coverage.png` - Boxplot across coverage levels

## Output Sections

### 1. Summary Statistics Table

Markdown table with bootstrap 95% confidence intervals:

| Category | Subcategory | N | Mean Error (%) | 95% CI | SD | Pass Rate (%) |
|----------|-------------|---|----------------|--------|----|--------------:|
| Overall | All tests | 270 | 2.74 | [2.53, 2.95] | 1.32 | 100.0 |
| Variant Type | SNP | 90 | 2.12 | [1.89, 2.35] | 0.98 | 100.0 |
| Variant Type | INS | 90 | 2.45 | [2.18, 2.72] | 1.15 | 100.0 |
| Variant Type | DEL | 90 | 2.53 | [2.26, 2.80] | 1.21 | 100.0 |

**Key Metrics**:
- **Mean Error**: Average percent error in recovering allelic ratios
- **95% CI**: Bootstrap confidence interval (10,000 iterations)
- **SD**: Standard deviation across replicates
- **Pass Rate**: Percentage of tests with error < 10%

### 2. Hypothesis Test: Unbiased Estimation

Tests if the algorithm has systematic bias:

- **H0**: Mean error = 0 (unbiased)
- **Ha**: Mean error ≠ 0 (biased)

**Example Output**:
```
Test: One-sample t-test
Statistic: 34.120
P-value: 1.23e-89
95% CI: [2.53, 2.95]%
Effect Size (Cohen's d): 2.07

Interpretation: REJECT H0 at α=0.05: Mean error significantly differs from 0
```

**What This Means**:
- Small positive bias (~2.7%) is statistically significant but well below 10% threshold
- Effect size is large (d=2.07), meaning the deviation from 0 is substantial in statistical terms
- However, the absolute error is still acceptably small for biological applications

### 3. Hypothesis Test: Performance Threshold

Tests if the algorithm meets our acceptance criteria:

- **H0**: Mean error ≥ 10% (algorithm FAILS)
- **Ha**: Mean error < 10% (algorithm PASSES)

**Example Output**:
```
Test: One-sample t-test (H0: μ ≥ 10%)
Statistic: -90.523
P-value: 3.45e-185
95% CI: [2.53, 2.95]%
Effect Size: 5.52 SDs below threshold

Interpretation: PASS: Mean error is significantly below 10% threshold
```

**What This Means**:
- Algorithm performs significantly better than our 10% threshold
- The mean error is 5.5 standard deviations below threshold (huge margin)
- p < 0.001 provides very strong evidence

### 4. Consistency and Reproducibility

**Metrics Reported**:

- **Coefficient of Variation (CV)**: Ratio of SD to mean
  - CV < 0.5: Low variability (excellent)
  - CV 0.5-1.0: Moderate variability (good)
  - CV > 1.0: High variability (concerning)

- **Intraclass Correlation (ICC)**: Measures reproducibility across replicates
  - ICC > 0.75: Excellent reproducibility
  - ICC 0.60-0.75: Good reproducibility
  - ICC < 0.60: Moderate reproducibility

**Example Output**:
```
Mean Error: 2.74%
Standard Deviation: 1.32%
Coefficient of Variation (CV): 0.48
Intraclass Correlation (ICC): 0.82

Interpretation: CV=0.48 indicates low variability.
                ICC=0.82 shows excellent reproducibility.
```

### 5. Sample Size Justification

Explains why N=10 replicates per configuration is sufficient:

**Example Output**:
```
Replicates per configuration: 10
Significance level (α): 0.05
Minimum Detectable Effect Size (MDES): 0.72
Achieved Power: 1.00

Justification: With N=10 replicates, we can detect effect sizes ≥0.72
               with 80% power at α=0.05. For a deterministic algorithm,
               this demonstrates consistency rather than statistical
               power estimation.
```

**Key Point**:
- We're testing a **deterministic algorithm**, not a statistical method
- Goal is to show **consistency**, not estimate power
- 10 replicates suffice to demonstrate algorithmic stability
- This follows WASP 2015 precedent (van de Geijn et al.)

### 6. Manuscript-Ready Text

The script generates copy-paste text for your manuscript:

#### Results Section (Example)

> Simulation validation (270 tests) demonstrated accurate recovery of planted allelic ratios with mean error of 2.7% (95% CI: [2.5, 3.0], SD=1.3%). One-sample t-test confirmed error was significantly below our 10% threshold (t=-90.5, p<0.001). Performance was consistent across variant types (SNP: 2.1% [1.9, 2.4], INS: 2.5% [2.2, 2.7], DEL: 2.5% [2.3, 2.8]). Coefficient of variation across 10 replicates per configuration was 0.48, indicating high algorithmic consistency.

#### Methods Section (Example)

> We validated WASP2 using metamorphic testing principles, verifying that: (1) planted allelic ratios are conserved within 10% error (conservation), (2) accuracy is consistent across variant types (symmetry), and (3) results are stable across random seeds (reproducibility). We tested 27 configurations with 10 replicates each (total n=270) to demonstrate consistency of the deterministic position mapping algorithm. This sample size provides adequate coverage of the parameter space for algorithmic validation, following the approach of van de Geijn et al. (2015). Statistical significance was assessed using one-sample t-tests with 95% confidence intervals calculated via bootstrap resampling (10,000 iterations).

## Command-Line Options

```bash
python analyze_simulation_results.py --help
```

### Options

- `results_file` (required): Path to `simulation_results.csv`
- `--output`, `-o`: Save markdown report to file
- `--plots`, `-p`: Directory for diagnostic plots
- `--alpha`: Significance level (default: 0.05)
- `--bootstrap`: Bootstrap iterations (default: 10,000)

### Examples

```bash
# Basic analysis (print to screen)
python analyze_simulation_results.py simulation_results.csv

# Save to file
python analyze_simulation_results.py simulation_results.csv -o report.md

# Generate plots
python analyze_simulation_results.py simulation_results.csv -p ./figures

# Stricter significance level
python analyze_simulation_results.py simulation_results.csv --alpha 0.01

# Full analysis with everything
python analyze_simulation_results.py simulation_results.csv \
    --output full_report.md \
    --plots ./figures \
    --alpha 0.05 \
    --bootstrap 20000
```

## Statistical Methods Explained

### Bootstrap Confidence Intervals

**What**: Non-parametric method for estimating uncertainty

**How**:
1. Resample data with replacement 10,000 times
2. Calculate statistic (mean, median, etc.) for each resample
3. Use percentiles of bootstrap distribution as CI bounds

**Why**:
- Doesn't assume normality
- Works for any statistic (mean, max, variance, etc.)
- More robust than parametric methods

**Implementation**:
```python
def bootstrap_confidence_interval(data, statistic=np.mean, n_bootstrap=10000):
    bootstrap_stats = []
    for i in range(n_bootstrap):
        sample = np.random.choice(data, size=len(data), replace=True)
        bootstrap_stats.append(statistic(sample))

    lower = np.percentile(bootstrap_stats, 2.5)
    upper = np.percentile(bootstrap_stats, 97.5)
    return (lower, upper)
```

### Hypothesis Testing Strategy

#### Test 1: Bias Detection
- **Purpose**: Detect systematic bias (should error be exactly 0?)
- **Test**: Two-sided t-test or Wilcoxon
- **Null**: Mean error = 0
- **Decision**: If p < 0.05, bias exists (but may be acceptable)

#### Test 2: Threshold Validation
- **Purpose**: Prove algorithm meets acceptance criteria
- **Test**: One-sided t-test
- **Null**: Mean error ≥ 10% (FAIL)
- **Decision**: If p < 0.05, algorithm PASSES

**Why Two Tests?**
- Test 1 detects *any* bias (even tiny)
- Test 2 proves bias is *acceptably small*
- Both can reject null (small bias that's still < 10%)

### Power Analysis for Sample Size

**Traditional Use**:
- Estimate sample size needed to detect effect with desired power
- Used when planning experiments

**Our Use**:
- Post-hoc justification of N=10
- Calculate minimum detectable effect size (MDES)
- Show 10 replicates are sufficient for consistency testing

**Key Formula**:
```
Power = P(reject H0 | H0 is false)
      = 1 - β (where β is Type II error rate)

MDES = effect size detectable with power=0.80 at α=0.05 for N replicates
```

**For N=10**:
- MDES ≈ 0.72 (Cohen's d)
- Can detect "medium" to "large" effects reliably
- Sufficient for deterministic algorithm validation

### Variance Metrics

#### Coefficient of Variation (CV)
```
CV = SD / Mean
```
- Normalized measure of variability
- CV < 0.5 is considered "low variability"
- Our CV ≈ 0.48 indicates consistent performance

#### Intraclass Correlation (ICC)
```
ICC = Between-group variance / (Between + Within variance)
```
- Measures consistency across replicates
- ICC > 0.75 indicates "excellent" reproducibility
- Our ICC ≈ 0.82 shows high reproducibility

## Interpreting Results

### ✅ Good Results

- Pass rate > 90%
- Mean error < 10%
- Threshold test: p < 0.05 (reject H0: μ ≥ 10%)
- CV < 0.5
- ICC > 0.75

### ⚠️ Warning Signs

- Pass rate 80-90%
- Mean error 5-10%
- High variance (CV > 0.7)
- Poor reproducibility (ICC < 0.60)

### ❌ Concerning Results

- Pass rate < 80%
- Mean error > 10%
- Threshold test: p > 0.05 (fail to reject H0)
- Very high variance (CV > 1.0)

## Addressing Reviewer Comments

### Common Questions and Answers

**Q: "Why only 10 replicates? Other studies use 1000+."**

**A**: We're testing a deterministic algorithm, not estimating statistical power. Our goal is to demonstrate consistency across different random seeds, not to characterize a distribution. Power analysis shows N=10 can detect effect sizes ≥0.72 with 80% power, which is sufficient for algorithmic validation. This follows the precedent of van de Geijn et al. (2015).

**Q: "What's the biological significance of a 2.7% error?"**

**A**: For allele-specific expression analysis, accurate allelic ratio estimation is critical. A 2.7% mean error translates to recovering a true 2:1 ratio as 1.95:1 or 2.05:1, which is biologically negligible and well within experimental variation. Our 10% threshold (recovering 2:1 as 1.8:1 to 2.2:1) is conservative.

**Q: "How do you know errors are normally distributed?"**

**A**: We test normality using Shapiro-Wilk test. If non-normal, we use non-parametric tests (Wilcoxon) instead of t-tests. Bootstrap CIs are non-parametric and don't assume normality.

**Q: "Did you correct for multiple testing?"**

**A**: We report overall error across all tests (primary analysis) and stratified analyses (secondary). For the primary hypothesis (overall mean < 10%), no correction is needed. For stratified analyses, we report descriptive statistics rather than inferential tests.

## Dependencies

```bash
pip install pandas numpy scipy matplotlib seaborn
```

Or use conda:
```bash
conda install pandas numpy scipy matplotlib seaborn
```

## Troubleshooting

### Error: "Missing required columns"

**Cause**: Results file doesn't have expected format

**Fix**: Ensure you're using output from `simulate_indel_ase_v2.py`

### Warning: "Filtered N tests with infinite error"

**Cause**: Some variants had zero ALT reads (division by zero)

**Impact**: These edge cases are excluded from analysis

**Fix**: Normal behavior - edge cases where no ALT reads aligned

### Error: "Too few samples for normality test"

**Cause**: Very small sample size (N < 3)

**Fix**: Use larger simulation tier (moderate or comprehensive)

## References

### Methods
- **Bootstrap**: Efron & Tibshirani (1993). *An Introduction to the Bootstrap*. Chapman & Hall/CRC.
- **Power Analysis**: Cohen (1988). *Statistical Power Analysis for the Behavioral Sciences*. Routledge.
- **Metamorphic Testing**: Chan et al. (2018). Testing Bioinformatics Programs. *IEEE Software*.

### WASP Validation
- van de Geijn et al. (2015). WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nature Methods*, 12(11), 1061-1063.

## Contact

For questions about statistical methods, see `CRITICAL_REVIEW_SIMULATION_APPROACH.md` for detailed justification of our approach.
