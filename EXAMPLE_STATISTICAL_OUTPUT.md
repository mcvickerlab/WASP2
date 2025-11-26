# Example Statistical Analysis Output

This document shows example output from `analyze_simulation_results.py` to demonstrate the type of publication-ready statistics you'll receive.

## Input Data

- **Simulation tier**: Moderate (270 tests)
- **Configurations**: 27 (3 variant types × 3 ratios × 3 coverage levels)
- **Replicates per config**: 10
- **Total tests**: 270

---

# Statistical Analysis of WASP2 Simulation Results

**Analysis Date**: 2025-11-25
**Total Tests**: 270
**Configurations**: 27
**Replicates per config**: 10

## Summary Statistics

| Category | Subcategory | N | Mean Error (%) | 95% CI | SD | Pass Rate (%) |
|----------|-------------|---|----------------|--------|----|--------------:|
| Overall | All tests | 270 | 2.74 | [2.53, 2.95] | 1.32 | 100.0 |
| Variant Type | SNP | 90 | 2.12 | [1.89, 2.35] | 0.98 | 100.0 |
| Variant Type | INS | 90 | 2.45 | [2.18, 2.72] | 1.15 | 100.0 |
| Variant Type | DEL | 90 | 2.53 | [2.26, 2.80] | 1.21 | 100.0 |
| Coverage | 20x | 90 | 3.21 | [2.85, 3.57] | 1.53 | 100.0 |
| Coverage | 50x | 90 | 2.74 | [2.45, 3.03] | 1.24 | 100.0 |
| Coverage | 100x | 90 | 2.27 | [2.01, 2.53] | 1.09 | 100.0 |
| Allelic Ratio | 1.0:1 | 90 | 2.65 | [2.38, 2.92] | 1.15 | 100.0 |
| Allelic Ratio | 2.0:1 | 90 | 2.71 | [2.43, 2.99] | 1.19 | 100.0 |
| Allelic Ratio | 4.0:1 | 90 | 2.86 | [2.56, 3.16] | 1.28 | 100.0 |

**Key Findings**:
- All 270 tests passed (100% pass rate)
- Overall mean error well below 10% threshold
- Performance consistent across variant types (SNP ≈ INS ≈ DEL)
- Error decreases with higher coverage (as expected)
- Error similar across allelic ratios (unbiased)

## Hypothesis Test: Unbiased Estimation

**Null Hypothesis (H0)**: Mean error = 0 (algorithm is unbiased)
**Alternative (Ha)**: Mean error ≠ 0 (algorithm has systematic bias)

**Test**: One-sample t-test
**Statistic**: 34.120
**P-value**: 1.23e-89
**95% CI**: [2.53, 2.95]%
**Effect Size (Cohen's d)**: 2.073

**Interpretation**: REJECT H0 at α=0.05: Mean error significantly differs from 0

**What This Means**:
The algorithm shows a small positive bias (~2.7%), which is statistically significant due to our large sample size (N=270) and low variance. However, this bias is **biologically negligible** and **well below our 10% acceptance threshold**. The effect size (d=2.07) indicates the deviation from zero is substantial in statistical terms, but the absolute magnitude (2.7%) is acceptable for practical applications.

**Biological Context**:
- A 2.7% bias means a true 2:1 ratio is recovered as ~1.95:1
- This is within typical biological/experimental variation
- Far better than naive read counting (which has ~20-30% reference bias)

## Hypothesis Test: Performance Threshold

**Null Hypothesis (H0)**: Mean error ≥ 10% (algorithm FAILS)
**Alternative (Ha)**: Mean error < 10% (algorithm PASSES)

**Test**: One-sample t-test (H0: μ ≥ 10%)
**Statistic**: -90.523
**P-value**: 3.45e-185
**95% CI for mean**: [2.53, 2.95]%
**Effect Size**: 5.52 SDs below threshold

**Interpretation**: PASS: Mean error is significantly below 10% threshold (p=3.45e-185)

**What This Means**:
- **STRONG EVIDENCE** that algorithm performs well
- Mean error is 5.5 standard deviations below our acceptance threshold
- P-value essentially zero (p < 10^-180)
- No reasonable doubt that algorithm meets validation criteria

**For Reviewers**:
This test provides rigorous statistical proof that WASP2's indel handling is accurate. The enormous test statistic (t=-90.5) and infinitesimal p-value indicate this is not a marginal result - the algorithm decisively passes validation.

## Consistency and Reproducibility

**Mean Error**: 2.74%
**Standard Deviation**: 1.32%
**Coefficient of Variation (CV)**: 0.48
**Intraclass Correlation (ICC)**: 0.82

**Interpretation**: CV=0.48 indicates low variability. ICC=0.82 shows excellent reproducibility.

**What This Means**:

1. **Low Variability (CV=0.48)**:
   - Normalized variation is less than 50% of the mean
   - Indicates consistent performance across different configurations
   - CV < 0.5 is considered "excellent" in most fields

2. **Excellent Reproducibility (ICC=0.82)**:
   - 82% of variance comes from true differences between configurations
   - Only 18% is random variation across replicates
   - ICC > 0.75 indicates "excellent" reproducibility
   - Shows algorithmic stability across random seeds

**Comparison to Benchmarks**:
- Biological assays: CV typically 0.5-1.5
- Genomic pipelines: CV typically 0.3-0.8
- Our results (CV=0.48): Within expected range, on the good side

## Sample Size Justification

**Replicates per configuration**: 10
**Significance level (α)**: 0.05
**Minimum Detectable Effect Size (MDES)**: 0.723
**Achieved Power**: 1.000

**Justification**: With N=10 replicates, we can detect effect sizes ≥0.72 with 80% power at α=0.05. For a deterministic algorithm, this demonstrates consistency rather than statistical power estimation.

**Detailed Explanation**:

### Why N=10 is Sufficient

1. **Deterministic Algorithm**:
   - WASP2's position mapping is deterministic (same input → same output)
   - Unlike statistical methods, we don't need 1000+ replicates to estimate distributions
   - Replicates test stability across different random seeds, not algorithm variability

2. **Effect Size Context**:
   - MDES of 0.72 corresponds to a "medium-to-large" effect (Cohen, 1988)
   - Our observed effect size (5.52) is 7.6× larger than MDES
   - Massively overpowered for our validation purposes

3. **Literature Precedent**:
   - WASP 2015 (van de Geijn et al.) used similar targeted simulation approach
   - BWA validation used 9 metamorphic relations (qualitative)
   - We test 27 configurations × 10 reps = more comprehensive than precedents

4. **Parameter Space Coverage**:
   - 3 variant types (SNP, INS, DEL) ✓
   - 3 coverage levels (20x, 50x, 100x) ✓
   - 3 allelic ratios (1:1, 2:1, 4:1) ✓
   - 10 replicates for consistency ✓
   - **Total: 270 tests** - excellent parameter space coverage

### What Would More Replicates Give Us?

- **N=10**: MDES = 0.72, CI width ≈ 0.42%
- **N=100**: MDES = 0.23, CI width ≈ 0.13%
- **N=1000**: MDES = 0.07, CI width ≈ 0.04%

With N=10, our 95% CI is [2.53, 2.95] (width = 0.42%).
Increasing to N=100 would narrow CI to ~[2.66, 2.79] (width = 0.13%).

**Diminishing Returns**:
- We already prove mean error < 10% with p < 10^-180
- Narrower CI doesn't change biological interpretation
- 10-fold increase in computational time for marginal gain

## Manuscript-Ready Text

### Results Section

Simulation validation (270 tests) demonstrated accurate recovery of planted allelic ratios with mean error of 2.7% (95% CI: [2.5, 3.0], SD=1.3%). One-sample t-test confirmed error was significantly below our 10% threshold (t=-90.5, p<0.001). Performance was consistent across variant types (SNP: 2.1% [1.9, 2.4], INS: 2.5% [2.2, 2.7], DEL: 2.5% [2.3, 2.8]) and coverage levels (20×: 3.2% [2.9, 3.6], 50×: 2.7% [2.5, 3.0], 100×: 2.3% [2.0, 2.5]). Coefficient of variation across 10 replicates per configuration was 0.48, indicating high algorithmic consistency.

### Methods Section

We validated WASP2 using metamorphic testing principles, verifying that: (1) planted allelic ratios are conserved within 10% error (conservation), (2) accuracy is consistent across variant types (symmetry), and (3) results are stable across random seeds (reproducibility). We tested 27 configurations (3 variant types × 3 allelic ratios × 3 coverage levels) with 10 replicates each (total n=270) to demonstrate consistency of the deterministic position mapping algorithm. This sample size provides adequate coverage of the parameter space for algorithmic validation, following the approach of van de Geijn et al. (2015). Statistical significance was assessed using one-sample t-tests with 95% confidence intervals calculated via bootstrap resampling (10,000 iterations).

## Diagnostic Plots

The following plots are generated with `--plots` option:

### 1. Error Distribution
- Shows histogram of all error values
- Red dashed line: 10% threshold
- Blue solid line: Mean error
- **Interpretation**: Tight distribution around 2.7%, all values < 10%

### 2. Q-Q Plot
- Tests normality assumption
- Points should fall on diagonal line if normal
- **Interpretation**: Slight deviation at tails, but acceptable for t-test

### 3. Error by Variant Type
- Boxplots comparing SNP, INS, DEL
- **Interpretation**: Similar median and IQR across types (symmetry)

### 4. Error by Coverage
- Boxplots comparing 20x, 50x, 100x
- **Interpretation**: Error decreases with coverage (expected)

## References

- van de Geijn et al. (2015). WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nature Methods*, 12(11), 1061-1063.
- Efron & Tibshirani (1993). An Introduction to the Bootstrap. *Chapman & Hall/CRC*.
- Cohen (1988). Statistical Power Analysis for the Behavioral Sciences. *Routledge*.

---

## How to Use This in Your Manuscript

### In Results Section

Copy the "Results Section" text above directly, or adapt as needed:

> Simulation validation (270 tests) demonstrated accurate recovery of planted allelic ratios...

### In Methods Section

Copy the "Methods Section" text above:

> We validated WASP2 using metamorphic testing principles...

### In Supplementary Materials

Include:
- **Supplementary Table S1**: Full summary statistics table (from above)
- **Supplementary Figure S1**: Error distribution plots
- **Supplementary Note S1**: Sample size justification (from above)

### Addressing Reviewer Comments

If reviewers ask:

**"Why only 10 replicates?"**
→ Point to Sample Size Justification section

**"Is the bias statistically significant?"**
→ Yes (p<0.001), but biologically negligible (2.7% vs 10% threshold)

**"How do you know it's reproducible?"**
→ ICC=0.82 indicates excellent reproducibility

**"What about multiple testing?"**
→ Primary hypothesis (overall error < 10%) tested once; stratified analyses are descriptive

## Summary: What You Can Claim

✅ **VALID CLAIMS**:
- "WASP2 accurately recovers planted allelic ratios (mean error 2.7%, 95% CI [2.5, 3.0])"
- "Performance significantly exceeds acceptance criteria (p < 0.001)"
- "Results are highly reproducible across replicates (ICC=0.82)"
- "Validation is consistent with metamorphic testing principles"
- "Sample size provides adequate parameter space coverage"

❌ **AVOID CLAIMING**:
- "WASP2 is perfectly unbiased" (we detect small 2.7% bias)
- "Results are normally distributed" (check Q-Q plot first)
- "Higher coverage doesn't affect accuracy" (it does - error decreases)

⚠️ **NUANCED CLAIMS**:
- "Small positive bias (2.7%) is statistically significant but biologically negligible"
- "N=10 replicates sufficient for deterministic algorithm validation (not power estimation)"

---

**Generated by**: `analyze_simulation_results.py`
**Based on**: CRITICAL_REVIEW_SIMULATION_APPROACH.md recommendations
