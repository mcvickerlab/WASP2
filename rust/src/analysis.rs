/// WASP2 Analysis Module - Beta-binomial Allelic Imbalance Detection
///
/// Rust implementation of the Python analysis stage (src/analysis/as_analysis.py)
/// Uses beta-binomial model to detect allelic imbalance in ASE data.
///
/// Performance target: 3-5x speedup over Python (2.7s → 0.5-0.9s)
use anyhow::{Context, Result};
use rayon::prelude::*;
use rv::dist::BetaBinomial;
use rv::traits::HasDensity;
use statrs::distribution::{ChiSquared, ContinuousCDF};
use std::collections::HashMap;

// ============================================================================
// Data Structures
// ============================================================================

/// Allele count data for a single variant
#[derive(Debug, Clone)]
pub struct VariantCounts {
    pub chrom: String,
    pub pos: u32,
    pub ref_count: u32,
    pub alt_count: u32,
    pub region: String,
}

/// Statistical results for a region
#[derive(Debug, Clone)]
pub struct ImbalanceResult {
    pub region: String,
    pub ref_count: u32,
    pub alt_count: u32,
    pub n: u32,
    pub snp_count: usize,
    pub null_ll: f64,  // Null model log-likelihood
    pub alt_ll: f64,   // Alternative model log-likelihood
    pub mu: f64,       // Estimated imbalance proportion
    pub lrt: f64,      // Likelihood ratio test statistic
    pub pval: f64,     // P-value
    pub fdr_pval: f64, // FDR-corrected p-value
}

/// Configuration for analysis
#[derive(Debug, Clone)]
pub struct AnalysisConfig {
    pub min_count: u32,
    pub pseudocount: u32,
    pub method: AnalysisMethod,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisMethod {
    Single, // Single dispersion parameter
    Linear, // Linear dispersion model
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            min_count: 10,
            pseudocount: 1,
            method: AnalysisMethod::Single,
        }
    }
}

// ============================================================================
// Rho Parameter Bounds (Issue #228)
// ============================================================================

/// Epsilon for clamping rho parameter to avoid division by zero.
/// Matches Python's RHO_EPSILON in as_analysis.py.
const RHO_EPSILON: f64 = 1e-10;

/// Clamp rho to safe range (epsilon, 1-epsilon) to prevent division by zero.
///
/// The beta-binomial parameterization uses alpha = mu * (1-rho) / rho, which
/// causes division by zero when rho=0 and produces zero alpha/beta when rho=1.
#[inline]
fn clamp_rho(rho: f64) -> f64 {
    rho.clamp(RHO_EPSILON, 1.0 - RHO_EPSILON)
}

// ============================================================================
// Core Statistical Functions
// ============================================================================

/// Calculate beta-binomial log-likelihood (negative for optimization)
///
/// Python equivalent: `opt_prob()` in as_analysis.py
///
/// # Arguments
/// * `prob` - Probability parameter (0 to 1)
/// * `rho` - Dispersion parameter (0 to 1), will be clamped to safe range
/// * `k` - Reference allele count
/// * `n` - Total count
///
/// # Returns
/// Negative log-likelihood value (for minimization)
pub fn opt_prob(prob: f64, rho: f64, k: u32, n: u32) -> Result<f64> {
    // Clamp rho to prevent division by zero (Issue #228)
    let rho = clamp_rho(rho);

    // Convert to alpha/beta parameters for beta-binomial
    let alpha = prob * (1.0 - rho) / rho;
    let beta = (1.0 - prob) * (1.0 - rho) / rho;

    // Create beta-binomial distribution (rv uses: n as u32, alpha, beta)
    let bb =
        BetaBinomial::new(n, alpha, beta).context("Failed to create beta-binomial distribution")?;

    // Return negative log-likelihood (rv uses reference for ln_f, k as u64)
    let log_pmf = bb.ln_f(&(k as u64));
    Ok(-log_pmf)
}

/// Calculate beta-binomial log-likelihood for array of counts
///
/// Python equivalent: Used in `single_model()` for null/alt likelihood
pub fn betabinom_logpmf_sum(
    ref_counts: &[u32],
    n_array: &[u32],
    alpha: f64,
    beta: f64,
) -> Result<f64> {
    let mut sum = 0.0;

    for (k, n) in ref_counts.iter().zip(n_array.iter()) {
        let bb = BetaBinomial::new(*n, alpha, beta)
            .context("Failed to create beta-binomial distribution")?;
        sum += bb.ln_f(&(*k as u64));
    }

    Ok(sum)
}

// ============================================================================
// Optimization Functions
// ============================================================================

/// Optimize dispersion parameter using Brent's method
///
/// Python equivalent: `minimize_scalar()` in scipy.optimize
fn optimize_dispersion(ref_counts: &[u32], n_array: &[u32]) -> Result<f64> {
    // Objective function: negative log-likelihood of null model (prob=0.5)
    let objective = |rho: f64| -> f64 {
        // Clamp rho to prevent division by zero (Issue #228)
        let rho = clamp_rho(rho);
        let alpha = 0.5 * (1.0 - rho) / rho;
        let beta = 0.5 * (1.0 - rho) / rho;

        match betabinom_logpmf_sum(ref_counts, n_array, alpha, beta) {
            Ok(ll) => -ll, // Return negative for minimization
            Err(_) => f64::INFINITY,
        }
    };

    // Use golden section search (simple but effective)
    let result = golden_section_search(objective, 0.001, 0.999, 1e-6)?;
    Ok(result)
}

/// Optimize probability parameter for alternative model
///
/// Python equivalent: `parse_opt()` calling `minimize_scalar(opt_prob, ...)`
fn optimize_prob(ref_counts: &[u32], n_array: &[u32], disp: f64) -> Result<(f64, f64)> {
    // For single SNP, optimize directly
    if ref_counts.len() == 1 {
        let objective = |prob: f64| -> f64 {
            match opt_prob(prob, disp, ref_counts[0], n_array[0]) {
                Ok(nll) => nll,
                Err(_) => f64::INFINITY,
            }
        };

        let mu = golden_section_search(objective, 0.0, 1.0, 1e-6)?;
        let alt_ll = -objective(mu);
        return Ok((alt_ll, mu));
    }

    // For multiple SNPs, sum log-likelihoods
    let objective = |prob: f64| -> f64 {
        let mut sum = 0.0;
        for (k, n) in ref_counts.iter().zip(n_array.iter()) {
            match opt_prob(prob, disp, *k, *n) {
                Ok(nll) => sum += nll,
                Err(_) => return f64::INFINITY,
            }
        }
        sum
    };

    let mu = golden_section_search(objective, 0.0, 1.0, 1e-6)?;
    let alt_ll = -objective(mu);
    Ok((alt_ll, mu))
}

/// Golden section search for 1D optimization
///
/// Simple but robust method for bounded scalar optimization.
/// Equivalent to scipy's minimize_scalar with method='bounded'
#[allow(unused_assignments)]
fn golden_section_search<F>(f: F, a: f64, mut b: f64, tol: f64) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    const PHI: f64 = 1.618033988749895; // Golden ratio
    let inv_phi = 1.0 / PHI;
    let inv_phi2 = 1.0 / (PHI * PHI);

    let mut a = a;
    let mut h = b - a;

    // Initial points
    let mut c = a + inv_phi2 * h;
    let mut d = a + inv_phi * h;
    let mut fc = f(c);
    let mut fd = f(d);

    // Iterate until convergence
    while h.abs() > tol {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            h = inv_phi * h;
            c = a + inv_phi2 * h;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            h = inv_phi * h;
            d = a + inv_phi * h;
            fd = f(d);
        }
    }

    Ok(if fc < fd { c } else { d })
}

// ============================================================================
// FDR Correction
// ============================================================================

/// Benjamini-Hochberg FDR correction
///
/// Python equivalent: `false_discovery_control(pvals, method="bh")`
pub fn fdr_correction(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    if n == 0 {
        return vec![];
    }

    // Create indexed p-values for sorting
    let mut indexed_pvals: Vec<(usize, f64)> = pvals.iter().copied().enumerate().collect();

    // Sort by p-value (ascending)
    indexed_pvals.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Calculate BH-adjusted p-values
    let mut adjusted = vec![0.0; n];
    let mut prev_adj = 1.0;

    for (rank, (idx, pval)) in indexed_pvals.iter().enumerate().rev() {
        let adj_pval = (pval * n as f64 / (rank + 1) as f64).min(prev_adj).min(1.0);
        adjusted[*idx] = adj_pval;
        prev_adj = adj_pval;
    }

    adjusted
}

// ============================================================================
// Main Analysis Functions
// ============================================================================

/// Single dispersion model analysis
///
/// Python equivalent: `single_model()` in as_analysis.py
pub fn single_model(variants: Vec<VariantCounts>) -> Result<Vec<ImbalanceResult>> {
    if variants.is_empty() {
        return Ok(vec![]);
    }

    // Extract ref_counts and N for all variants
    let ref_counts: Vec<u32> = variants.iter().map(|v| v.ref_count).collect();
    let n_array: Vec<u32> = variants.iter().map(|v| v.ref_count + v.alt_count).collect();

    // Step 1: Optimize global dispersion parameter
    eprintln!("Optimizing dispersion parameter...");
    let disp = optimize_dispersion(&ref_counts, &n_array)?;
    eprintln!("  Dispersion: {:.6}", disp);

    // Step 2: Group by region
    let mut region_map: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, variant) in variants.iter().enumerate() {
        region_map
            .entry(variant.region.clone())
            .or_default()
            .push(i);
    }

    eprintln!(
        "Optimizing imbalance likelihood for {} regions...",
        region_map.len()
    );

    // Step 3: Calculate null and alternative likelihoods per region (parallel)
    // Clamp disp before calculating null_param (Issue #228)
    let disp = clamp_rho(disp);
    let null_param = 0.5 * (1.0 - disp) / disp;

    let results: Result<Vec<_>> = region_map
        .par_iter()
        .map(|(region, indices)| -> Result<ImbalanceResult> {
            // Extract counts for this region
            let region_ref: Vec<u32> = indices.iter().map(|&i| ref_counts[i]).collect();
            let region_n: Vec<u32> = indices.iter().map(|&i| n_array[i]).collect();

            // Null model: prob = 0.5 (no imbalance)
            let null_ll = betabinom_logpmf_sum(&region_ref, &region_n, null_param, null_param)?;

            // Alternative model: optimize prob
            let (alt_ll, mu) = optimize_prob(&region_ref, &region_n, disp)?;

            // Likelihood ratio test
            let lrt = -2.0 * (null_ll - alt_ll);

            // P-value from chi-squared distribution (df=1)
            let chi2 = ChiSquared::new(1.0).context("Failed to create chi-squared distribution")?;
            let pval = 1.0 - chi2.cdf(lrt);

            // Sum counts for this region
            let total_ref: u32 = region_ref.iter().sum();
            let total_alt: u32 = indices.iter().map(|&i| variants[i].alt_count).sum();
            let total_n = total_ref + total_alt;

            Ok(ImbalanceResult {
                region: region.clone(),
                ref_count: total_ref,
                alt_count: total_alt,
                n: total_n,
                snp_count: indices.len(),
                null_ll,
                alt_ll,
                mu,
                lrt,
                pval,
                fdr_pval: 0.0, // Will be filled later
            })
        })
        .collect();

    let mut results = results?;

    // Step 4: FDR correction
    let pvals: Vec<f64> = results.iter().map(|r| r.pval).collect();
    let fdr_pvals = fdr_correction(&pvals);

    for (result, fdr_pval) in results.iter_mut().zip(fdr_pvals.iter()) {
        result.fdr_pval = *fdr_pval;
    }

    Ok(results)
}

/// Main entry point for allelic imbalance analysis
///
/// Python equivalent: `get_imbalance()` in as_analysis.py
pub fn analyze_imbalance(
    variants: Vec<VariantCounts>,
    config: &AnalysisConfig,
) -> Result<Vec<ImbalanceResult>> {
    // Apply filters and pseudocounts
    let filtered: Vec<VariantCounts> = variants
        .into_iter()
        .map(|mut v| {
            v.ref_count += config.pseudocount;
            v.alt_count += config.pseudocount;
            v
        })
        .filter(|v| {
            let n = v.ref_count + v.alt_count;
            n >= config.min_count + (2 * config.pseudocount)
        })
        .collect();

    eprintln!("Processing {} variants after filtering", filtered.len());

    // Run analysis based on method
    let mut results = match config.method {
        AnalysisMethod::Single => single_model(filtered)?,
        AnalysisMethod::Linear => {
            return Err(anyhow::anyhow!("Linear model not yet implemented"));
        }
    };

    // Remove pseudocounts from results
    for result in results.iter_mut() {
        if result.ref_count < config.pseudocount
            || result.alt_count < config.pseudocount
            || result.n < 2 * config.pseudocount
        {
            eprintln!(
                "[WARN] Counts smaller than pseudocount for region {}: ref={}, alt={}, n={}, pc={}",
                result.region, result.ref_count, result.alt_count, result.n, config.pseudocount
            );
        }
        result.ref_count = result.ref_count.saturating_sub(config.pseudocount);
        result.alt_count = result.alt_count.saturating_sub(config.pseudocount);
        result.n = result.n.saturating_sub(2 * config.pseudocount);
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_opt_prob() {
        // Test beta-binomial likelihood calculation
        let result = opt_prob(0.5, 0.1, 10, 20).unwrap();
        assert!(result.is_finite());
        assert!(result > 0.0); // Negative log-likelihood should be positive
    }

    #[test]
    fn test_opt_prob_rho_boundary_zero() {
        // Issue #228: rho=0 should not cause division by zero
        let result = opt_prob(0.5, 0.0, 10, 20).unwrap();
        assert!(result.is_finite(), "rho=0 should produce finite result after clamping");
        assert!(!result.is_nan(), "rho=0 should not produce NaN");
    }

    #[test]
    fn test_opt_prob_rho_boundary_one() {
        // Issue #228: rho=1 should not produce zero alpha/beta
        let result = opt_prob(0.5, 1.0, 10, 20).unwrap();
        assert!(result.is_finite(), "rho=1 should produce finite result after clamping");
        assert!(!result.is_nan(), "rho=1 should not produce NaN");
    }

    #[test]
    fn test_opt_prob_rho_near_boundaries() {
        // Test values very close to boundaries
        for rho in [1e-15, 1e-12, 1e-10, 0.999999999, 0.9999999999999] {
            let result = opt_prob(0.5, rho, 10, 20).unwrap();
            assert!(result.is_finite(), "rho={} should produce finite result", rho);
        }
    }

    #[test]
    fn test_clamp_rho() {
        // Test clamping function directly
        assert!((clamp_rho(0.0) - RHO_EPSILON).abs() < 1e-15);
        assert!((clamp_rho(1.0) - (1.0 - RHO_EPSILON)).abs() < 1e-15);
        assert!((clamp_rho(0.5) - 0.5).abs() < 1e-15);
        assert!((clamp_rho(-1.0) - RHO_EPSILON).abs() < 1e-15);
        assert!((clamp_rho(2.0) - (1.0 - RHO_EPSILON)).abs() < 1e-15);
    }

    #[test]
    fn test_fdr_correction() {
        let pvals = vec![0.01, 0.05, 0.1, 0.5];
        let fdr = fdr_correction(&pvals);

        // FDR-adjusted p-values should be >= original
        for (orig, adj) in pvals.iter().zip(fdr.iter()) {
            assert!(adj >= orig);
        }
    }

    #[test]
    fn test_golden_section() {
        // Test optimization on simple quadratic
        let f = |x: f64| (x - 0.7).powi(2);
        let min = golden_section_search(f, 0.0, 1.0, 1e-6).unwrap();
        assert!((min - 0.7).abs() < 1e-5);
    }
}
