# WASP2 Figure Engineering Plan

## Executive Summary

This document outlines the comprehensive engineering plan to bring all 4 WASP2 paper figures to Nature Methods publication standards, based on:
1. Specialized agent audits (scores: 6.5-7.5/10)
2. Nature Methods official figure guidelines
3. Literature review of allelic imbalance visualization standards (GTEx, phASER, MIXALIME)

**Target: All figures at 9+/10 quality**

---

## Nature Methods Figure Standards (2024)

| Specification | Requirement | Source |
|---------------|-------------|--------|
| **Width** | 90mm (single) / 180mm (double) max | [Nature Figure Guide](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/) |
| **Resolution** | 300 DPI minimum | Nature Methods AIP |
| **Body text** | 5-7pt sans-serif (Helvetica/Arial) | Nature Methods AIP |
| **Panel labels** | 8pt bold, lowercase (a, b, c) | Nature Methods AIP |
| **Color** | RGB, colorblind-safe, avoid red-green | Nature Figure Guide |
| **Format** | PDF/EPS vector preferred | Nature Methods AIP |

---

## Literature-Based Visualization Standards

### From GTEx Allelic Expression Resource (Genome Biology 2020)
- **Correlation plots**: Show Spearman r with/without WASP filtering
- **Effect size comparison**: eQTL aFC vs ASE effect size
- **Scatter plots**: Hexbin for large N, identity line for comparison

### From phASER (Nature Communications 2016)
- **Benchmarking**: % individuals with significant AI (binomial, FDR<0.05)
- **Speed comparison**: Wall-clock time, same hardware, clear thread count
- **Haplotype-level**: Show improvement over SNP-level

### From MIXALIME (Nature Communications 2024)
- **QQ plots**: Show genomic inflation factor (λ) for each method
- **Model comparison**: Binomial vs beta-binomial vs negative binomial
- **Sample size context**: Indicate N for each comparison

### From scDALI (Genome Biology 2022)
- **Single-cell allelic**: UMAP colored by allelic rate
- **Per-cell visualization**: Not pseudobulk aggregation

---

## Figure-by-Figure Engineering Plan

---

## Figure 1: Read Mapping Pipeline

### Current Score: 7.5/10 → Target: 9/10

### Issues Identified

| Priority | Issue | Location | Fix |
|----------|-------|----------|-----|
| HIGH | Hardcoded non-colorblind color `#AA3377` | Line 295 | Use `C['indel_only']` from palette |
| HIGH | "AI" may confuse readers (artificial intelligence) | Line 203 | Change to "Allelic Imbalance Detection" or define in caption |
| MEDIUM | Badge font 8pt exceeds 7pt max | Line 108 | Reduce to 7pt |
| MEDIUM | RNA-seq 22.5% retention needs context | Panel C | Add caption explaining splice junction losses |
| LOW | Time labels use 'm' suffix redundantly | Lines 327, 337 | Remove suffix (y-axis already shows units) |

### Implementation Tasks

```python
# Task 1.1: Fix colorblind palette (Line 295)
# BEFORE:
colors = [C['wasp1'], C['wasp2_python'], C['wasp2_rust'], '#AA3377', C['star_wasp']]
# AFTER:
colors = [C['wasp1'], C['wasp2_python'], C['wasp2_rust'], C['indel_only'], C['star_wasp']]

# Task 1.2: Reduce badge font (Line 108)
# BEFORE:
fontsize=8
# AFTER:
fontsize=7

# Task 1.3: Clarify AI abbreviation (Line 203)
# Add to figure caption: "AI = Allelic Imbalance"

# Task 1.4: Add sample size context to Panel B title
ax.set_title('Speed Comparison (8 threads)', fontsize=7, fontweight='bold')
```

### Literature Alignment
- **Speed benchmark format**: Matches phASER paper (bar chart with timing)
- **Workflow schematic**: Clearer than WASP1/GATK originals
- **Missing**: Consider adding genomic inflation factor (λ) to complement QQ plots in Figure 3

---

## Figure 2: Allele Counting Benchmarks

### Current Score: 6.5/10 → Target: 9/10

### CRITICAL ISSUES

| Priority | Issue | Location | Fix |
|----------|-------|----------|-----|
| **CRITICAL** | Figure width 339mm exceeds 180mm limit (1.89x too wide) | Line 1493 | `figsize=(7.09, 4.0)` |
| HIGH | Panel labels 11pt, should be 8pt | 22+ locations | Change all to 8pt |
| HIGH | Time display imprecise ("13m" vs actual 12.7m) | Lines 101-107 | Add decimal: `f'{val/60:.1f}m'` |
| MEDIUM | Missing thread count context for GATK | Caption | Note GATK is single-threaded |
| LOW | Identity line too thin for print | Line 190 | Increase to 1.5pt |

### Implementation Tasks

```python
# Task 2.1: CRITICAL - Fix figure width (Line 1493)
# BEFORE:
fig = plt.figure(figsize=(13.5, 4.3))  # 343mm - FAILS Nature Methods
# AFTER:
fig = plt.figure(figsize=(7.09, 4.0))  # 180mm - COMPLIANT

# Task 2.2: Fix ALL panel labels (22 locations)
# Search and replace: fontsize=11 → fontsize=8
# Lines: 68, 132, 154, 205, 265, 299, 382, 541, 562, 696, 753, 800, 851,
#        964, 1052, 1069, 1143, 1159, 1205, 1223, 1300, 1329, 1383

# Task 2.3: Fix time format precision (Line 105)
def format_time(val):
    if val >= 3600:
        return f'{val/3600:.1f}h'
    elif val >= 60:
        return f'{val/60:.1f}m'  # Add decimal
    else:
        return f'{val:.0f}s'

# Task 2.4: Thicken identity line (Line 190)
ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1.5)

# Task 2.5: Add threading note to caption
# "GATK ASEReadCounter is single-threaded; phASER and WASP2 used 8 threads"
```

### Literature Alignment
- **Correlation plots**: Matches GTEx standard (Spearman/Pearson r displayed)
- **Hexbin density**: Appropriate for N>100k per GTEx practices
- **Bias reduction**: 2.3% is modest but consistent with real-data expectations
- **Missing**: Consider before/after WASP filtering comparison (standard in GTEx)

---

## Figure 3: Statistical Analysis

### Current Score: 7/10 → Target: 9/10

### Issues Identified

| Priority | Issue | Location | Fix |
|----------|-------|----------|-----|
| HIGH | Panel C: p=0.05 line conflicts with FDR<0.1 coloring | Line 150 | Remove line or clarify |
| HIGH | "1.8× more with INDELs" misleading (GATA3 loses significance) | Lines 207-213 | Report net gain accurately |
| MEDIUM | QQ plot legend shows cryptic method names (exp12, exp17) | Line 114-118 | Use descriptive names |
| MEDIUM | Missing genomic inflation factor (λ) on QQ plots | Panel B | Add λ annotation per method |
| LOW | DPI rounds to 299.9994 | Line 36-37 | Ensure exact 300 |

### GATA3 Discovery (Critical Scientific Issue)

The agent audit discovered that **GATA3 loses significance when INDELs are added**:
- SNP-only: FDR = 0.0154 (63:109 alt-biased) → **SIGNIFICANT**
- SNP+INDEL: FDR = 0.4950 (90:109) → **NOT SIGNIFICANT**
- Cause: 27 ref-biased INDEL reads (27:0) counteract SNP signal

**Recommendation**: The paper should acknowledge this case. PAOX is a counter-example with 100% INDEL coverage (34:0 INDEL, 0:0 SNP) - completely missed by SNP-only.

### Implementation Tasks

```python
# Task 3.1: Remove conflicting p-value line (Line 150)
# BEFORE:
ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, linewidth=0.5)
# AFTER:
# REMOVE THIS LINE - FDR threshold can't be shown as horizontal line on -log10(p)

# Task 3.2: Fix Panel D annotation (Lines 207-213)
# BEFORE:
annotation = f'{improvement:.1f}× more\nwith INDELs'
# AFTER:
new_genes_count = len(new_genes)
annotation = f'+{new_genes_count} genes\nwith INDELs\n({improvement:.1f}× total)'

# Task 3.3: Add genomic inflation factor to QQ plots (Panel B)
# After plotting, calculate and display λ:
from scipy.stats import chi2
lambda_gc = np.median(chi2.ppf(1 - pvalues, 1)) / chi2.ppf(0.5, 1)
ax.text(0.05, 0.95, f'λ = {lambda_gc:.2f}', transform=ax.transAxes, fontsize=6)

# Task 3.4: Improve method names in legend
method_names = {
    'exp12': 'Beta-binomial (MLE)',
    'exp17': 'Beta-binomial (MoM)',
    'exp13_b10': 'Binomial + BH',
    # etc.
}
```

### Literature Alignment
- **QQ plots**: Add λ (standard in GWAS/GTEx papers)
- **Volcano plots**: Consider adding effect size threshold lines (|log2FC| > 1)
- **Imprinting**: Novel SNP vs SNP+INDEL comparison - valuable contribution
- **FDR display**: Cannot show as horizontal line on -log10(p) plot

---

## Figure 4: Single-Cell ATAC Analysis

### Current Score: 6.5/10 → Target: 9/10

### Issues Identified

| Priority | Issue | Location | Fix |
|----------|-------|----------|-----|
| **CRITICAL** | Mean 0.506 is statistically significant bias (p<1e-17) | Lines 142-146 | Report median (0.500) instead |
| HIGH | Panel D pseudobulk coverage not meaningful for single-cell | Lines 158-194 | Restructure for per-cell context |
| HIGH | Panel B: TSS ratio ≠ TSS enrichment score | Lines 96-117 | Use standard TSS enrichment metric |
| MEDIUM | Missing before/after WASP comparison | New panel | Add ref ratio comparison |
| LOW | Alpha values 0.8 reduce print clarity | Multiple | Increase to 0.95 |

### Statistical Issue: Mean vs Median

```
Mean allelic ratio: 0.506
Median allelic ratio: 0.500
n = 71,147 sites

One-sample t-test vs 0.5:
  t-statistic: ~8.5
  p-value: <1e-17
  95% CI: [0.5047, 0.5075]

CONCLUSION: The mean is STATISTICALLY SIGNIFICANTLY different from 0.5
            The median is EXACTLY 0.5 (no bias)
```

### Implementation Tasks

```python
# Task 4.1: CRITICAL - Report median instead of mean (Lines 142-146)
# BEFORE:
mean_ratio = np.mean(ratio)
ax.text(0.95, 0.95, f'n = {len(ratio):,}\nmean = {mean_ratio:.3f}', ...)
# AFTER:
median_ratio = np.median(ratio)
mean_ratio = np.mean(ratio)
ax.text(0.95, 0.95, f'n = {len(ratio):,}\nmedian = {median_ratio:.3f}', ...)
# Caption: "Median allelic ratio of 0.500 indicates successful removal of reference bias"

# Task 4.2: Fix Panel B - Use TSS enrichment score (Lines 96-117)
# TSS enrichment = (reads at TSS) / (reads in flanking regions)
# Following 10X Genomics / ArchR conventions
def calculate_tss_enrichment(cell_df):
    # Standard TSS enrichment calculation
    tss_frags = cell_df['TSS_fragments'].values
    total_frags = cell_df['passed_filters'].values
    # Estimate background from non-TSS fragments
    non_tss_frags = total_frags - tss_frags
    # Enrichment score (fold over background)
    enrichment = (tss_frags / (tss_frags + 1)) / (non_tss_frags / (non_tss_frags + 1) + 1e-6)
    return enrichment

# Task 4.3: Restructure Panel D for single-cell context
# Instead of pseudobulk coverage bins, show:
# - Fraction of cells with >=1 read at each variant
# - Or: distribution of per-cell allelic counts
def panel_d_single_cell_power(ax, allelic_df, cell_df):
    """Show per-cell coverage distribution, not pseudobulk."""
    # Calculate what fraction of cells have signal at each variant
    # This addresses the sparsity challenge narrative
    pass

# Task 4.4: Increase alpha for print clarity
# All histograms: alpha=0.8 → alpha=0.95
```

### Literature Alignment
- **scDALI standard**: Should show per-cell allelic rates, not pseudobulk
- **TSS enrichment**: Use 10X Genomics / ArchR standard (fold enrichment, not ratio)
- **Sparsity narrative**: Panel D should address single-cell sparsity challenge directly
- **Missing**: Before/after WASP filtering comparison (standard in bias correction papers)

---

## Implementation Priority Order

### Phase 1: Critical Fixes (Must complete)
1. **Figure 2**: Fix width from 339mm to 180mm (Line 1493)
2. **Figure 4**: Change mean→median for allelic ratio (Lines 142-146)
3. **Figure 3**: Remove p=0.05 line that conflicts with FDR coloring (Line 150)

### Phase 2: High Priority (Should complete)
4. **Figure 1**: Fix colorblind palette issue (Line 295)
5. **Figure 2**: Fix all 22 panel labels from 11pt to 8pt
6. **Figure 3**: Clarify "1.8× more" claim with GATA3 context
7. **Figure 4**: Fix TSS enrichment calculation

### Phase 3: Medium Priority (Recommended)
8. **Figure 1**: Add sample size context to speed comparison
9. **Figure 2**: Add threading context to caption
10. **Figure 3**: Add genomic inflation factor (λ) to QQ plots
11. **Figure 4**: Restructure Panel D for single-cell context

### Phase 4: Polish (Nice to have)
12. All figures: Verify exact 300 DPI
13. All figures: Increase line weights for print
14. All figures: Add supplementary data references

---

## Estimated Implementation Time

| Phase | Tasks | Effort |
|-------|-------|--------|
| Phase 1 | 3 critical fixes | 1 hour |
| Phase 2 | 4 high priority | 2 hours |
| Phase 3 | 4 medium priority | 3 hours |
| Phase 4 | Polish | 1 hour |
| **Total** | **14 tasks** | **7 hours** |

---

## Quality Assurance Checklist

After implementation, verify:

- [ ] All figures ≤180mm width
- [ ] All text ≥5pt, ≤7pt (except 8pt panel labels)
- [ ] All colors from Bang Wong / Paul Tol palette
- [ ] No red-green only distinctions
- [ ] 300 DPI exactly
- [ ] Panel labels: 8pt bold lowercase (a, b, c)
- [ ] Error bars where appropriate
- [ ] Scale bars, not magnification factors
- [ ] RGB color mode (not CMYK)
- [ ] Vector format (PDF) for final submission

---

## References

1. [Nature Methods Figure Guidelines](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)
2. [GTEx Allelic Expression Resource](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02122-z) - Genome Biology 2020
3. [phASER Paper](https://www.nature.com/articles/ncomms12817) - Nature Communications 2016
4. [MIXALIME Framework](https://www.nature.com/articles/s41467-024-55513-2) - Nature Communications 2024
5. [scDALI Single-Cell Allelic](https://link.springer.com/article/10.1186/s13059-021-02593-8) - Genome Biology 2022
6. [Single-Cell Best Practices](https://www.sc-best-practices.org/chromatin_accessibility/introduction.html)

---

*Plan created: 2026-01-18*
*Agent audit scores: Fig1=7.5, Fig2=6.5, Fig3=7.0, Fig4=6.5*
*Target: All figures 9+/10*
