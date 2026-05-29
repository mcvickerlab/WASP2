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
use std::collections::{HashMap, HashSet};

// ============================================================================
// Data Structures
// ============================================================================

/// Allele count data for a single variant
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VariantCounts {
    pub chrom: String,
    pub pos: u32,
    pub ref_count: u32,
    pub alt_count: u32,
    pub region: String,
    /// Donor/sample identifier; required only by the per_donor model.
    pub sample: Option<String>,
    /// Genotype phase relative to the reference allele: 0 = ref|alt, 1 = alt|ref.
    /// `None` when genotypes are absent or unphased. Required by phased analysis.
    pub gt: Option<u8>,
}

/// Deduplication key for matching Python's `df[keep_cols].drop_duplicates()`.
/// Used to remove exact duplicate observations before analysis.
#[derive(Hash, Eq, PartialEq)]
struct DedupeKey {
    chrom: String,
    pos: u32,
    region: String,
    ref_count: u32,
    alt_count: u32,
    n: u32,
    gt: Option<u8>, // Include GT when phased=true
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
    pub phased: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisMethod {
    Single,   // Single global dispersion parameter
    Linear,   // Linear dispersion model: rho = expit(d1 + N*d2)
    PerDonor, // Per-donor dispersion: rho fit separately per sample
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            min_count: 10,
            pseudocount: 1,
            method: AnalysisMethod::Single,
            phased: false,
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
///
/// # Note on Silent Clamping
///
/// This function does not log warnings when clamping occurs because it is called
/// in tight optimization loops where logging would impact performance. Extreme
/// rho values (at boundaries) may indicate data quality issues. The Python
/// implementation (`as_analysis.py`) supports optional warning via `warn=True`.
#[inline]
fn clamp_rho(rho: f64) -> f64 {
    rho.clamp(RHO_EPSILON, 1.0 - RHO_EPSILON)
}

/// Deduplicate variants to match Python's `df[keep_cols].drop_duplicates()`.
///
/// Python keeps columns: chrom, pos, ref_count, alt_count, N, region (+ GT if phased)
/// and removes exact duplicate rows. This ensures parity with Python's behavior.
fn deduplicate_variants(variants: Vec<VariantCounts>, phased: bool) -> Vec<VariantCounts> {
    let mut seen: HashSet<DedupeKey> = HashSet::new();
    let mut result = Vec::with_capacity(variants.len());

    for v in variants {
        let key = DedupeKey {
            chrom: v.chrom.clone(),
            pos: v.pos,
            region: v.region.clone(),
            ref_count: v.ref_count,
            alt_count: v.alt_count,
            n: v.ref_count + v.alt_count,
            gt: if phased { v.gt } else { None },
        };
        if seen.insert(key) {
            result.push(v);
        }
    }
    result
}

/// Sigmoid function (inverse logit): expit(x) = 1 / (1 + e^-x)
///
/// Numerically stable implementation that avoids overflow for large |x|.
#[inline]
fn expit(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
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

/// Optimize linear dispersion parameters using Nelder-Mead
///
/// Python equivalent: `minimize(opt_linear, x0=(0, 0), method="Nelder-Mead")`
/// Fits ρ = expit(d1 + N × d2) to the data.
fn optimize_dispersion_linear(ref_counts: &[u32], n_array: &[u32]) -> Result<(f64, f64)> {
    let objective = |params: [f64; 2]| -> f64 {
        let (d1, d2) = (params[0], params[1]);
        let mut neg_ll = 0.0;

        for (&ref_count, &n) in ref_counts.iter().zip(n_array.iter()) {
            // Clamp to prevent overflow in expit
            let exp_in = (d1 + n as f64 * d2).clamp(-10.0, 10.0);
            let rho = clamp_rho(expit(exp_in));
            let alpha = 0.5 * (1.0 - rho) / rho;

            // Beta-binomial log-pmf with α = β (null: p=0.5)
            if let Ok(bb) = rv::dist::BetaBinomial::new(n, alpha, alpha) {
                use rv::traits::HasDensity;
                neg_ll -= bb.ln_f(&(ref_count as u64));
            } else {
                return f64::INFINITY;
            }
        }

        neg_ll
    };

    // Multi-start Nelder-Mead. A single start at (0,0) can stall on the flat plateau
    // created by the exp_in clamp (when |d1 + N*d2| > 10 for all N, rho is constant, so
    // the surface has no descent direction and the simplex collapses at an unconverged,
    // worse point). Run NM from a neutral grid of starts spanning the plausible (d1,d2)
    // range and keep the best (lowest neg-LL). Deterministic; starts are NOT seeded near
    // any expected answer. Fixes the case where the simple simplex returned e.g. (1.5,-2.0).
    let starts: [[f64; 2]; 8] = [
        [0.0, 0.0],
        [-3.0, 0.0],
        [-6.0, 0.0],
        [3.0, 0.0],
        [0.0, -0.2],
        [-3.0, -0.2],
        [-6.0, -0.2],
        [3.0, -0.2],
    ];
    let mut best: Option<(f64, (f64, f64))> = None;
    for s in starts {
        if let Ok((d1, d2)) = nelder_mead_2d(&objective, s, 1e-6, 1000) {
            let val = objective([d1, d2]);
            if val.is_finite() && best.map_or(true, |(bv, _)| val < bv) {
                best = Some((val, (d1, d2)));
            }
        }
    }
    best.map(|(_, p)| p)
        .ok_or_else(|| anyhow::anyhow!("linear dispersion optimization failed from all starts"))
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

/// Optimize probability for phased data using known genotype phase.
///
/// Python equivalent: `parse_opt()` phased branch calling
/// `minimize_scalar(opt_phased_new, ...)` in as_analysis.py (L213-256, L118-140).
///
/// The phase array is first normalized so the first SNP is the reference
/// (`if gt[0] > 0 { gt = 1 - gt }`), matching Python's
/// `if gt_array[0] > 0: gt_array = 1 - gt_array`. The objective then sums
/// `opt_prob(|prob - gt_i|, disp, ref_i, n_i)` across the region and is
/// minimized over `prob` in (0, 1).
///
/// # Arguments
/// * `ref_counts` - Reference allele counts for the region
/// * `n_array` - Total counts for the region
/// * `gt_array` - Genotype phase per variant (0 or 1)
/// * `disp` - Dispersion parameter (scalar; per-region as in single/linear alt path)
///
/// # Returns
/// Tuple of (alt_ll, mu) - alternative-model log-likelihood and imbalance proportion
fn optimize_prob_phased(
    ref_counts: &[u32],
    n_array: &[u32],
    gt_array: &[u8],
    disp: f64,
) -> Result<(f64, f64)> {
    // Normalize phase so the first SNP is with respect to ref (matches Python).
    let flip = gt_array.first().is_some_and(|&g| g > 0);
    let gt_norm: Vec<u8> = gt_array
        .iter()
        .map(|&g| if flip { 1 - g } else { g })
        .collect();

    // objective(prob) = Σ opt_prob(|prob - gt_i|, disp, ref_i, n_i)
    let objective = |prob: f64| -> f64 {
        let mut sum = 0.0;
        for ((&k, &n), &g) in ref_counts.iter().zip(n_array.iter()).zip(gt_norm.iter()) {
            let phased_prob = (prob - g as f64).abs();
            match opt_prob(phased_prob, disp, k, n) {
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

/// Optimize probability for unphased data using dynamic programming.
///
/// Python equivalent: `opt_unphased_dp()` in as_analysis.py
///
/// This marginalizes over unknown phase using dynamic programming. At each
/// position after the first, we consider both phase possibilities and weight
/// them equally (0.5 each), accumulating likelihood across the region.
///
/// # Arguments
/// * `prob` - Probability parameter (0 to 1)
/// * `disp` - Dispersion parameter(s). If slice length > 1, first element is for
///   first position and remaining elements are for subsequent positions.
/// * `first_ref` - Reference count for first position
/// * `first_n` - Total count for first position
/// * `phase_ref` - Reference counts for subsequent positions
/// * `phase_n` - Total counts for subsequent positions
///
/// # Returns
/// Negative log-likelihood value (for minimization)
pub fn optimize_prob_unphased_dp(
    prob: f64,
    disp: &[f64],
    first_ref: u32,
    first_n: u32,
    phase_ref: &[u32],
    phase_n: &[u32],
) -> Result<f64> {
    // Split per-variant dispersion for first vs subsequent positions (linear model)
    let (first_disp, phase_disp): (f64, &[f64]) = if disp.len() > 1 {
        (disp[0], &disp[1..])
    } else {
        // Scalar dispersion: use same value for all positions
        (disp[0], disp)
    };

    // Get likelihood of first position (negative log-likelihood)
    let first_ll = opt_prob(prob, first_disp, first_ref, first_n)?;

    // If no subsequent positions, just return first_ll
    if phase_ref.is_empty() {
        return Ok(first_ll);
    }

    // Compute likelihoods for subsequent positions under both phase hypotheses
    // phase1_like: pmf with prob
    // phase2_like: pmf with (1 - prob)
    let mut prev_like: f64 = 1.0;

    for (i, (&ref_count, &n)) in phase_ref.iter().zip(phase_n.iter()).enumerate() {
        // Get dispersion for this position
        let rho = if phase_disp.len() > 1 {
            clamp_rho(phase_disp[i.min(phase_disp.len() - 1)])
        } else {
            clamp_rho(phase_disp[0])
        };

        // Phase 1: use prob
        let alpha1 = prob * (1.0 - rho) / rho;
        let beta1 = (1.0 - prob) * (1.0 - rho) / rho;
        let bb1 = BetaBinomial::new(n, alpha1, beta1)
            .context("Failed to create beta-binomial for phase1")?;
        let p1 = bb1.f(&(ref_count as u64)); // pmf (not log)

        // Phase 2: use (1 - prob)
        let alpha2 = (1.0 - prob) * (1.0 - rho) / rho;
        let beta2 = prob * (1.0 - rho) / rho;
        let bb2 = BetaBinomial::new(n, alpha2, beta2)
            .context("Failed to create beta-binomial for phase2")?;
        let p2 = bb2.f(&(ref_count as u64)); // pmf (not log)

        // DP update: marginalize over both phases with equal weight
        let p1_combined = prev_like * p1;
        let p2_combined = prev_like * p2;
        prev_like = 0.5 * p1_combined + 0.5 * p2_combined;
    }

    // Return first_ll + (-ln(prev_like))
    // Handle edge case where prev_like might be 0 or very small
    if prev_like <= 0.0 || !prev_like.is_finite() {
        return Ok(f64::INFINITY);
    }

    Ok(first_ll + (-prev_like.ln()))
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

/// Nelder-Mead simplex optimization for 2D functions
///
/// Used for fitting linear dispersion parameters (d1, d2).
/// Equivalent to scipy's minimize(..., method="Nelder-Mead")
fn nelder_mead_2d<F>(f: F, x0: [f64; 2], tol: f64, max_iter: usize) -> Result<(f64, f64)>
where
    F: Fn([f64; 2]) -> f64,
{
    // Standard Nelder-Mead coefficients
    let alpha = 1.0; // reflection
    let gamma = 2.0; // expansion
    let rho = 0.5; // contraction
    let sigma = 0.5; // shrinkage

    // Initialize simplex: equilateral triangle around x0
    let step = 1.0;
    let mut simplex = [
        (x0, f(x0)),
        ([x0[0] + step, x0[1]], f([x0[0] + step, x0[1]])),
        ([x0[0], x0[1] + step], f([x0[0], x0[1] + step])),
    ];

    for _iter in 0..max_iter {
        // Sort by function value (ascending)
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let (x_best, f_best) = simplex[0];
        let (x_worst, f_worst) = simplex[2];
        let (x_second, f_second) = simplex[1];

        // Check convergence: simplex size
        let size = ((x_worst[0] - x_best[0]).powi(2) + (x_worst[1] - x_best[1]).powi(2)).sqrt();
        if size < tol {
            return Ok((x_best[0], x_best[1]));
        }

        // Centroid of all points except worst
        let centroid = [
            (x_best[0] + x_second[0]) / 2.0,
            (x_best[1] + x_second[1]) / 2.0,
        ];

        // Reflection
        let x_r = [
            centroid[0] + alpha * (centroid[0] - x_worst[0]),
            centroid[1] + alpha * (centroid[1] - x_worst[1]),
        ];
        let f_r = f(x_r);

        if f_r < f_second && f_r >= f_best {
            // Accept reflection
            simplex[2] = (x_r, f_r);
            continue;
        }

        if f_r < f_best {
            // Try expansion
            let x_e = [
                centroid[0] + gamma * (x_r[0] - centroid[0]),
                centroid[1] + gamma * (x_r[1] - centroid[1]),
            ];
            let f_e = f(x_e);

            if f_e < f_r {
                simplex[2] = (x_e, f_e);
            } else {
                simplex[2] = (x_r, f_r);
            }
            continue;
        }

        // f_r >= f_second: try contraction
        let x_c = if f_r < f_worst {
            // Outside contraction
            [
                centroid[0] + rho * (x_r[0] - centroid[0]),
                centroid[1] + rho * (x_r[1] - centroid[1]),
            ]
        } else {
            // Inside contraction
            [
                centroid[0] + rho * (x_worst[0] - centroid[0]),
                centroid[1] + rho * (x_worst[1] - centroid[1]),
            ]
        };
        let f_c = f(x_c);

        if f_c < f_worst {
            simplex[2] = (x_c, f_c);
            continue;
        }

        // Shrink toward best point
        for i in 1..3 {
            let (x_i, _) = simplex[i];
            let x_new = [
                x_best[0] + sigma * (x_i[0] - x_best[0]),
                x_best[1] + sigma * (x_i[1] - x_best[1]),
            ];
            simplex[i] = (x_new, f(x_new));
        }
    }

    // Return best point after max iterations
    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    Ok((simplex[0].0[0], simplex[0].0[1]))
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
pub fn single_model(variants: Vec<VariantCounts>, phased: bool) -> Result<Vec<ImbalanceResult>> {
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

            // Alternative model: optimize prob (phase-aware routing)
            let all_gt = indices.iter().all(|&i| variants[i].gt.is_some());
            let (alt_ll, mu) = if phased && indices.len() > 1 && all_gt {
                // Phased multi-SNP region: use known genotype phase.
                let region_gt: Vec<u8> = indices
                    .iter()
                    .map(|&i| variants[i].gt.expect("gt present (checked by all_gt)"))
                    .collect();
                optimize_prob_phased(&region_ref, &region_n, &region_gt, disp)?
            } else if !phased && indices.len() > 1 {
                // Unphased multi-SNP region: marginalize over phase via DP,
                // minimizing opt_unphased_dp over prob in (0, 1)
                // (first = index 0, phase = remaining), matching Python's
                // parse_opt() unphased branch (minimize_scalar(opt_unphased_dp,...)).
                let disp_slice = [disp];
                let first_ref = region_ref[0];
                let first_n = region_n[0];
                let phase_ref = &region_ref[1..];
                let phase_n = &region_n[1..];
                let objective = |prob: f64| -> f64 {
                    match optimize_prob_unphased_dp(
                        prob,
                        &disp_slice,
                        first_ref,
                        first_n,
                        phase_ref,
                        phase_n,
                    ) {
                        Ok(nll) => nll,
                        Err(_) => f64::INFINITY,
                    }
                };
                let mu = golden_section_search(objective, 0.0, 1.0, 1e-6)?;
                let alt_ll = -objective(mu);
                (alt_ll, mu)
            } else {
                optimize_prob(&region_ref, &region_n, disp)?
            };

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

    // Sort by region name for deterministic output across runs
    results.sort_by(|a, b| a.region.cmp(&b.region));

    // Step 4: FDR correction
    let pvals: Vec<f64> = results.iter().map(|r| r.pval).collect();
    let fdr_pvals = fdr_correction(&pvals);

    for (result, fdr_pval) in results.iter_mut().zip(fdr_pvals.iter()) {
        result.fdr_pval = *fdr_pval;
    }

    Ok(results)
}

/// Linear dispersion model analysis
///
/// Python equivalent: `linear_model()` in as_analysis.py
/// Uses ρ = expit(d1 + N × d2) where d1, d2 are globally optimized parameters.
pub fn linear_model(variants: Vec<VariantCounts>, phased: bool) -> Result<Vec<ImbalanceResult>> {
    if variants.is_empty() {
        return Ok(vec![]);
    }

    // Extract ref_counts and N for all variants
    let ref_counts: Vec<u32> = variants.iter().map(|v| v.ref_count).collect();
    let n_array: Vec<u32> = variants.iter().map(|v| v.ref_count + v.alt_count).collect();

    // Step 1: Optimize linear dispersion parameters (d1, d2)
    eprintln!("Optimizing linear dispersion parameters...");
    let (d1, d2) = optimize_dispersion_linear(&ref_counts, &n_array)?;
    eprintln!("  d1 = {:.6}, d2 = {:.6}", d1, d2);

    // Step 2: Compute per-variant rho array
    let rho_array: Vec<f64> = n_array
        .iter()
        .map(|&n| {
            let exp_in = (d1 + n as f64 * d2).clamp(-10.0, 10.0);
            clamp_rho(expit(exp_in))
        })
        .collect();

    // Step 3: Group by region
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

    // Step 4: Calculate null and alternative likelihoods per region (parallel)
    let results: Result<Vec<_>> = region_map
        .par_iter()
        .map(|(region, indices)| -> Result<ImbalanceResult> {
            // Extract counts and rho for this region
            let region_ref: Vec<u32> = indices.iter().map(|&i| ref_counts[i]).collect();
            let region_n: Vec<u32> = indices.iter().map(|&i| n_array[i]).collect();
            let region_rho: Vec<f64> = indices.iter().map(|&i| rho_array[i]).collect();

            // Null model: prob = 0.5, using per-variant rho
            let mut null_ll = 0.0;
            for ((&ref_count, &n), &rho) in region_ref
                .iter()
                .zip(region_n.iter())
                .zip(region_rho.iter())
            {
                let alpha = 0.5 * (1.0 - rho) / rho;
                let bb = rv::dist::BetaBinomial::new(n, alpha, alpha)
                    .context("Failed to create beta-binomial for null")?;
                use rv::traits::HasDensity;
                null_ll += bb.ln_f(&(ref_count as u64));
            }

            // Alternative model: optimize prob using per-variant rho
            // For simplicity, use average rho for optimization (matches Python behavior)
            let avg_rho = region_rho.iter().sum::<f64>() / region_rho.len() as f64;
            // Phase-aware routing (disp = avg_rho scalar, matching the existing
            // single/linear alt path).
            let all_gt = indices.iter().all(|&i| variants[i].gt.is_some());
            let (alt_ll, mu) = if phased && indices.len() > 1 && all_gt {
                // Phased multi-SNP region: use known genotype phase.
                let region_gt: Vec<u8> = indices
                    .iter()
                    .map(|&i| variants[i].gt.expect("gt present (checked by all_gt)"))
                    .collect();
                optimize_prob_phased(&region_ref, &region_n, &region_gt, avg_rho)?
            } else if !phased && indices.len() > 1 {
                // Unphased multi-SNP region: marginalize over phase via DP,
                // minimizing opt_unphased_dp over prob in (0, 1)
                // (first = index 0, phase = remaining), matching Python's
                // parse_opt() unphased branch (minimize_scalar(opt_unphased_dp,...)).
                let disp_slice = [avg_rho];
                let first_ref = region_ref[0];
                let first_n = region_n[0];
                let phase_ref = &region_ref[1..];
                let phase_n = &region_n[1..];
                let objective = |prob: f64| -> f64 {
                    match optimize_prob_unphased_dp(
                        prob,
                        &disp_slice,
                        first_ref,
                        first_n,
                        phase_ref,
                        phase_n,
                    ) {
                        Ok(nll) => nll,
                        Err(_) => f64::INFINITY,
                    }
                };
                let mu = golden_section_search(objective, 0.0, 1.0, 1e-6)?;
                let alt_ll = -objective(mu);
                (alt_ll, mu)
            } else {
                optimize_prob(&region_ref, &region_n, avg_rho)?
            };

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
                fdr_pval: 0.0,
            })
        })
        .collect();

    let mut results = results?;

    // Step 5: FDR correction
    let pvals: Vec<f64> = results.iter().map(|r| r.pval).collect();
    let fdr_pvals = fdr_correction(&pvals);

    for (result, fdr_pval) in results.iter_mut().zip(fdr_pvals.iter()) {
        result.fdr_pval = *fdr_pval;
    }

    Ok(results)
}

/// Fit per-donor dispersion (rho) via bounded MLE, one value per `sample`.
///
/// Python equivalent: `compute_per_donor_rho()` — each donor's rho is the
/// symmetric beta-binomial MLE over all that donor's observations, bounded to
/// (1e-6, 1-1e-6). No extra floor (matches locked donor_rho.tsv: min ~7e-6).
fn compute_donor_rho(variants: &[VariantCounts]) -> Result<HashMap<String, f64>> {
    // Group observation indices by donor
    let mut donor_map: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, v) in variants.iter().enumerate() {
        let donor = v.sample.clone().ok_or_else(|| {
            anyhow::anyhow!("per_donor model requires a 'sample' column in the counts TSV")
        })?;
        donor_map.entry(donor).or_default().push(i);
    }

    // Fit rho per donor in parallel (golden-section on symmetric beta-binomial)
    donor_map
        .par_iter()
        .map(|(donor, idx)| -> Result<(String, f64)> {
            let dref: Vec<u32> = idx.iter().map(|&i| variants[i].ref_count).collect();
            let dn: Vec<u32> = idx
                .iter()
                .map(|&i| variants[i].ref_count + variants[i].alt_count)
                .collect();
            let objective = |rho: f64| -> f64 {
                let rho = clamp_rho(rho);
                let a = 0.5 * (1.0 - rho) / rho;
                match betabinom_logpmf_sum(&dref, &dn, a, a) {
                    Ok(ll) => -ll,
                    Err(_) => f64::INFINITY,
                }
            };
            // Bounds (1e-6, 1-1e-6) match scipy minimize_scalar(method="bounded")
            let rho = golden_section_search(objective, 1e-6, 1.0 - 1e-6, 1e-6)?;
            Ok((donor.clone(), rho))
        })
        .collect()
}

/// Per-donor dispersion model analysis.
///
/// Python equivalent: WASP2 `per_donor` extension — fit rho per `sample`, then
/// run the standard per-region beta-binomial LRT using each observation's
/// donor-specific rho (per-observation rho in BOTH null and alternative).
pub fn per_donor_model(variants: Vec<VariantCounts>) -> Result<Vec<ImbalanceResult>> {
    if variants.is_empty() {
        return Ok(vec![]);
    }

    // Step 1: fit per-donor dispersion
    eprintln!("Fitting per-donor dispersion...");
    let donor_rho = compute_donor_rho(&variants)?;
    eprintln!("  Fit rho for {} donors", donor_rho.len());
    // Emit fitted rho (for validation against locked donor_rho.tsv)
    let mut donors_sorted: Vec<(&String, &f64)> = donor_rho.iter().collect();
    donors_sorted.sort_by(|a, b| a.0.cmp(b.0));
    for (donor, rho) in donors_sorted {
        eprintln!("DONOR_RHO\t{}\t{:.8}", donor, rho);
    }

    // Step 2: per-observation rho = each obs's donor rho
    let ref_counts: Vec<u32> = variants.iter().map(|v| v.ref_count).collect();
    let n_array: Vec<u32> = variants.iter().map(|v| v.ref_count + v.alt_count).collect();
    let rho_array: Vec<f64> = variants
        .iter()
        .map(|v| {
            let donor = v
                .sample
                .as_ref()
                .expect("sample present (checked in compute_donor_rho)");
            clamp_rho(
                *donor_rho
                    .get(donor)
                    .expect("donor rho fit for every sample"),
            )
        })
        .collect();

    // Step 3: group by region
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

    // Step 4: per-region null/alt LL + LRT, using per-observation rho
    let results: Result<Vec<_>> = region_map
        .par_iter()
        .map(|(region, indices)| -> Result<ImbalanceResult> {
            let region_ref: Vec<u32> = indices.iter().map(|&i| ref_counts[i]).collect();
            let region_n: Vec<u32> = indices.iter().map(|&i| n_array[i]).collect();
            let region_rho: Vec<f64> = indices.iter().map(|&i| rho_array[i]).collect();

            use rv::traits::HasDensity;

            // Null model: prob = 0.5, per-observation rho
            let mut null_ll = 0.0;
            for ((&ref_count, &n), &rho) in region_ref
                .iter()
                .zip(region_n.iter())
                .zip(region_rho.iter())
            {
                let alpha = 0.5 * (1.0 - rho) / rho;
                let bb = rv::dist::BetaBinomial::new(n, alpha, alpha)
                    .context("Failed to create beta-binomial for null")?;
                null_ll += bb.ln_f(&(ref_count as u64));
            }

            // Alternative model: optimize one pi, per-observation rho
            // (matches Python _alt_ll: a=pi(1-rho)/rho, b=(1-pi)(1-rho)/rho)
            let alt_obj = |pi: f64| -> f64 {
                let pi = pi.clamp(1e-6, 1.0 - 1e-6);
                let mut nll = 0.0;
                for ((&ref_count, &n), &rho) in region_ref
                    .iter()
                    .zip(region_n.iter())
                    .zip(region_rho.iter())
                {
                    let a = pi * (1.0 - rho) / rho;
                    let b = (1.0 - pi) * (1.0 - rho) / rho;
                    match rv::dist::BetaBinomial::new(n, a, b) {
                        Ok(bb) => nll -= bb.ln_f(&(ref_count as u64)),
                        Err(_) => return f64::INFINITY,
                    }
                }
                nll
            };
            let mu = golden_section_search(alt_obj, 1e-6, 1.0 - 1e-6, 1e-6)?;
            let alt_ll = -alt_obj(mu);

            // Likelihood ratio test (clamped at 0 like Python)
            let lrt = (-2.0 * (null_ll - alt_ll)).max(0.0);
            let chi2 = ChiSquared::new(1.0).context("Failed to create chi-squared distribution")?;
            let pval = 1.0 - chi2.cdf(lrt);

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
                fdr_pval: 0.0,
            })
        })
        .collect();

    let mut results = results?;

    // Deterministic output order
    results.sort_by(|a, b| a.region.cmp(&b.region));

    // FDR correction (BH)
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

    // Deduplicate variants for Python parity (single and linear models only).
    // Python's get_imbalance() does: df[keep_cols].drop_duplicates()
    // Per-donor model skips dedup as it's a Rust-only feature.
    let deduped = if config.method != AnalysisMethod::PerDonor {
        let before = filtered.len();
        let after = deduplicate_variants(filtered, config.phased);
        if after.len() < before {
            eprintln!(
                "Deduplicated {} → {} variants ({} duplicates removed)",
                before,
                after.len(),
                before - after.len()
            );
        }
        after
    } else {
        filtered
    };

    // Run analysis based on method
    let mut results = match config.method {
        AnalysisMethod::Single => single_model(deduped, config.phased)?,
        AnalysisMethod::Linear => linear_model(deduped, config.phased)?,
        AnalysisMethod::PerDonor => per_donor_model(deduped)?,
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
        assert!(
            result.is_finite(),
            "rho=0 should produce finite result after clamping"
        );
        assert!(!result.is_nan(), "rho=0 should not produce NaN");
    }

    #[test]
    fn test_opt_prob_rho_boundary_one() {
        // Issue #228: rho=1 should not produce zero alpha/beta
        let result = opt_prob(0.5, 1.0, 10, 20).unwrap();
        assert!(
            result.is_finite(),
            "rho=1 should produce finite result after clamping"
        );
        assert!(!result.is_nan(), "rho=1 should not produce NaN");
    }

    #[test]
    fn test_opt_prob_rho_near_boundaries() {
        // Test values very close to boundaries
        for rho in [1e-15, 1e-12, 1e-10, 0.999999999, 0.9999999999999] {
            let result = opt_prob(0.5, rho, 10, 20).unwrap();
            assert!(
                result.is_finite(),
                "rho={} should produce finite result",
                rho
            );
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

    #[test]
    fn test_optimize_prob_unphased_dp_single_position() {
        // Test with only first position (no subsequent positions)
        let disp = vec![0.1];
        let result = optimize_prob_unphased_dp(0.5, &disp, 10, 20, &[], &[]).unwrap();
        // Should just return first_ll from opt_prob
        let expected = opt_prob(0.5, 0.1, 10, 20).unwrap();
        assert!(
            (result - expected).abs() < 1e-10,
            "Single position: result={}, expected={}",
            result,
            expected
        );
    }

    #[test]
    fn test_optimize_prob_unphased_dp_multiple_positions() {
        // Test with multiple positions
        let disp = vec![0.1];
        let phase_ref = vec![8, 12];
        let phase_n = vec![20, 25];
        let result = optimize_prob_unphased_dp(0.5, &disp, 10, 20, &phase_ref, &phase_n).unwrap();
        // Result should be finite and positive
        assert!(result.is_finite(), "Result should be finite");
        assert!(result > 0.0, "Negative log-likelihood should be positive");
    }

    #[test]
    fn test_optimize_prob_unphased_dp_per_variant_disp() {
        // Test with per-variant dispersion (linear model)
        let disp = vec![0.1, 0.15, 0.2]; // first pos: 0.1, subsequent: [0.15, 0.2]
        let phase_ref = vec![8, 12];
        let phase_n = vec![20, 25];
        let result = optimize_prob_unphased_dp(0.5, &disp, 10, 20, &phase_ref, &phase_n).unwrap();
        assert!(
            result.is_finite(),
            "Result should be finite with per-variant disp"
        );
        assert!(result > 0.0, "Negative log-likelihood should be positive");
    }

    #[test]
    fn test_optimize_prob_unphased_dp_asymmetric_prob() {
        // Test with asymmetric probability (imbalance)
        let disp = vec![0.1];
        let phase_ref = vec![15, 18]; // More ref alleles
        let phase_n = vec![20, 25];
        let result_balanced =
            optimize_prob_unphased_dp(0.5, &disp, 10, 20, &phase_ref, &phase_n).unwrap();
        let result_imbalanced =
            optimize_prob_unphased_dp(0.7, &disp, 14, 20, &phase_ref, &phase_n).unwrap();
        // Both should be finite
        assert!(result_balanced.is_finite());
        assert!(result_imbalanced.is_finite());
    }
}
