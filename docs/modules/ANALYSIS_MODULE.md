# Analysis Module - Technical Overview

**Phase 1.3 Deep Dive**
**Total LOC:** 2,779 (48.1% of codebase)
**Files:** 10 Python files
**Purpose:** Statistical analysis of allelic imbalance using beta-binomial models

---

## Module Structure

### Core Statistical Engine
- **as_analysis.py** (674 LOC) - Beta-binomial optimization and likelihood calculations
- **as_analysis_sc.py** (258 LOC) - Single-cell statistical methods

### Orchestrators
- **run_analysis.py** (204 LOC) - Bulk analysis pipeline
- **run_analysis_sc.py** (266 LOC) - Single-cell analysis pipeline
- **run_compare_ai.py** (77 LOC) - Differential analysis orchestrator

### Comparison Engine
- **compare_ai.py** (516 LOC) - Compare allelic imbalance between cell types/conditions

### Data Processing
- **filter_data.py** (124 LOC) - VCF/GTF filtering, BAM processing
- **count_alleles.py** (122 LOC) - Bulk allele counting (analysis version)
- **count_alleles_sc.py** (186 LOC) - Single-cell allele counting (analysis version)

### CLI
- **__main__.py** (352 LOC) - 3 commands: find_imbalance, find_imbalance_sc, compare_imbalance

---

## Statistical Methodology

### Beta-Binomial Model

**Null Hypothesis (H0):** No allelic imbalance (μ = 0.5)
```python
# Dispersion optimization
opt_disp = lambda rho, ref, n: -np.sum(
    betabinom.logpmf(ref, n, (0.5 * (1 - rho) / rho), (0.5 * (1 - rho) / rho)))
disp = minimize_scalar(opt_disp, bounds=(0,1))
```

**Alternative Hypothesis (H1):** Allelic imbalance exists (μ ≠ 0.5)
```python
# Optimize probability parameter
mu = minimize_scalar(opt_prob, args=(disp, ref, n), bounds=(0, 1))["x"]
```

**Likelihood Ratio Test:**
```python
lrt = -2 * (null_ll - alt_ll)
pval = chi2.sf(lrt, 1)  # Chi-squared with 1 degree of freedom
```

### Two Dispersion Models

1. **Single Model** - One dispersion parameter for entire dataset
2. **Linear Model** - Weighted dispersion by total read count N

### Phased vs Unphased

**Phased (with genotype info):**
```python
def opt_phased_new(prob, disp, ref_data, n_data, gt_data):
    phased_ll = opt_prob(np.abs(prob - gt_data), disp, ref_data, n_data)
    return np.sum(phased_ll)
```

**Unphased (dynamic programming):**
```python
def opt_unphased_dp(prob, disp, first_ref, first_n, phase_ref, phase_n):
    first_ll = opt_prob(prob, disp, first_ref[0], first_n[0])
    # Accumulate likelihoods across phase combinations
    for p1, p2 in zip(phase1_like, phase2_like):
        prev_like = (0.5 * p1_combined) + (0.5 * p2_combined)
    return first_ll + -np.log(prev_like)
```

---

## Pipeline Flow

### Bulk Analysis (find_imbalance)
1. Load count TSV
2. Add pseudocounts and filter by min_count
3. Estimate dispersion parameter (ρ)
4. Optimize μ per region/SNP
5. Calculate LRT and p-values
6. FDR correction (Benjamini-Hochberg)
7. Write results TSV

### Single-Cell Analysis (find_imbalance_sc)
1. Load h5ad AnnData file
2. Parse barcode mappings → cell type groups
3. Filter heterozygous genotypes
4. QC: Remove outlier SNPs/regions (z-score cutoff)
5. Global dispersion estimation across all data
6. **Per cell type:**
   - Aggregate counts across cells in group
   - Filter regions by min_count
   - Optimize μ per region
   - LRT and FDR correction
7. Write separate TSV per cell type

### Differential Analysis (compare_imbalance)
1. Load h5ad file
2. Process counts per cell type
3. **For each pair of cell types:**
   - Find shared SNPs in regions
   - H0: Same μ for both groups (combined optimization)
   - H1: Different μ (separate optimizations)
   - LRT comparing null vs alternative
4. Write TSV per comparison

---

## Key Classes

### WaspAnalysisData (Bulk)
```python
class WaspAnalysisData:
    count_file, min_count, pseudocount, phased, model,
    out_file, region_col, groupby
```

### WaspAnalysisSC (Single-Cell)
```python
class WaspAnalysisSC:
    adata_file, bc_map, min_count, pseudocount, phased,
    sample, groups, model, out_file, z_cutoff
```

---

## Critical Design Patterns

### 1. Nested Optimization
```python
# Null: Combined likelihood
null_res = minimize_scalar(opt_combined_imbalance, args=(disp, func1, func2))

# Alt: Separate likelihoods
alt_res1 = minimize_scalar(like_func1, args=(disp, ...))
alt_res2 = minimize_scalar(like_func2, args=(disp, ...))
```

### 2. Dynamic Likelihood Selection
```python
def get_imbalance_func(ref_count, n_count, phase_array=None):
    if len(ref_count) == 1:
        return opt_prob, (ref_count[0], n_count[0])
    elif phase_array is None:
        return opt_unphased_dp, (first_ref, first_n, phase_ref, phase_n)
    else:
        return opt_phased_new, (ref_count, n_count, phase_array)
```

### 3. AnnData Processing Pipeline
```python
process_adata_inputs() → adata_count_qc() → get_imbalance_sc() → write results
```

---

## Data Structures

### Input: Count TSV (Bulk)
```
chr  pos  ref  alt  [GT]  [region]  [parent]  ref_count  alt_count  other_count
```

### Input: h5ad AnnData (Single-Cell)
```python
adata.obs        # SNP metadata (chr, pos, ref, alt, GT)
adata.var        # Barcodes + group assignments
adata.layers     # "ref" and "alt" sparse count matrices
adata.uns        # "feature" mapping (SNP → region), "samples" list
```

### Output: Results TSV
```
region  num_snps  mu  null_ll  alt_ll  pval  fdr_pval
```

### Output: Comparison TSV
```
region  num_snps  combined_mu  mu1  mu2  null_ll  alt_ll  pval  fdr_pval
```

---

## Module Complexity Metrics

| File | LOC | Functions | Classes | Complexity |
|------|-----|-----------|---------|------------|
| as_analysis.py | 674 | 9 | 0 | Very High |
| compare_ai.py | 516 | 6 | 0 | High |
| __main__.py | 352 | 3 | 0 | Medium |
| run_analysis_sc.py | 266 | 2 | 1 | Medium |
| as_analysis_sc.py | 258 | 3 | 0 | High |
| run_analysis.py | 204 | 2 | 1 | Low |
| count_alleles_sc.py | 186 | 5 | 0 | Medium |
| filter_data.py | 124 | 7 | 0 | Low |
| count_alleles.py | 122 | 3 | 0 | Low |
| run_compare_ai.py | 77 | 1 | 0 | Low |

---

## Dependencies by Function

| Dependency | Usage |
|------------|-------|
| numpy | Array operations, optimization inputs |
| pandas | DataFrame manipulation, CSV I/O |
| scipy.stats | betabinom, chi2, false_discovery_control |
| scipy.optimize | minimize_scalar for μ and ρ optimization |
| anndata | h5ad file I/O for single-cell data |
| pysam | BAM/VCF reading, pileup operations |
| pybedtools | BED/VCF/GTF intersections |

---

## See Also

- **ANALYSIS_ISSUES.md** - Technical debt and code quality issues
- **COUNTING_MODULE.md** - Upstream count generation
- **ARCHITECTURE.md** - Overall system design
