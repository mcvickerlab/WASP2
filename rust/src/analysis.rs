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

/// Results and fitted nuisance parameters from one independent analysis table.
#[derive(Debug, Clone)]
pub struct AnalysisOutput {
    pub results: Vec<ImbalanceResult>,
    pub rho: Option<f64>,
    pub linear_params: Option<(f64, f64)>,
    pub linear_depth_center: Option<f64>,
    pub linear_depth_scale: Option<f64>,
    pub n_observations: usize,
    pub effective_phased: bool,
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
    Single, // Single global dispersion parameter
    Linear, // Linear dispersion model: rho = expit(d1 + d2*((N-center)/scale))
}

/// Nuisance parameters fitted under the balanced beta-binomial null.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DispersionParameters {
    Single {
        rho: f64,
    },
    Linear {
        d1: f64,
        d2: f64,
        depth_center: f64,
        depth_scale: f64,
    },
}

impl DispersionParameters {
    fn validate_for(self, method: AnalysisMethod) -> Result<Self> {
        match (method, self) {
            (AnalysisMethod::Single, Self::Single { rho }) => {
                if !rho.is_finite() || !(0.0..=1.0).contains(&rho) {
                    anyhow::bail!("fixed scalar rho must be finite and in [0, 1]");
                }
            }
            (
                AnalysisMethod::Linear,
                Self::Linear {
                    d1,
                    d2,
                    depth_center,
                    depth_scale,
                },
            ) => {
                if !d1.is_finite() || !d2.is_finite() {
                    anyhow::bail!("fixed linear dispersion parameters must be finite");
                }
                if !depth_center.is_finite() {
                    anyhow::bail!("fixed linear depth center must be finite");
                }
                if !depth_scale.is_finite() || depth_scale <= 0.0 {
                    anyhow::bail!("fixed linear depth scale must be finite and positive");
                }
            }
            (AnalysisMethod::Single, Self::Linear { .. }) => {
                anyhow::bail!("method 'single' requires a fixed scalar rho");
            }
            (AnalysisMethod::Linear, Self::Single { .. }) => {
                anyhow::bail!(
                    "method 'linear' requires fixed linear_d1, linear_d2, \
                     linear_depth_center, and linear_depth_scale"
                );
            }
        }
        Ok(self)
    }

    #[inline]
    fn rho_at_depth(self, n: u32) -> f64 {
        match self {
            Self::Single { rho } => clamp_rho(rho),
            Self::Linear {
                d1,
                d2,
                depth_center,
                depth_scale,
            } => {
                let standardized_depth = (n as f64 - depth_center) / depth_scale;
                clamp_rho(expit(d1 + standardized_depth * d2))
            }
        }
    }

    pub fn rho(self) -> Option<f64> {
        match self {
            Self::Single { rho } => Some(rho),
            Self::Linear { .. } => None,
        }
    }

    pub fn linear_params(self) -> Option<(f64, f64)> {
        match self {
            Self::Single { .. } => None,
            Self::Linear { d1, d2, .. } => Some((d1, d2)),
        }
    }

    pub fn linear_depth_standardization(self) -> Option<(f64, f64)> {
        match self {
            Self::Single { .. } => None,
            Self::Linear {
                depth_center,
                depth_scale,
                ..
            } => Some((depth_center, depth_scale)),
        }
    }
}

/// Result of nuisance-only fitting on one independent table.
#[derive(Debug, Clone, Copy)]
pub struct DispersionFit {
    pub parameters: DispersionParameters,
    pub n_observations: usize,
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

fn effective_phased(variants: &[VariantCounts], requested_phased: bool) -> bool {
    if !requested_phased {
        return false;
    }

    let mut saw_ref_alt = false;
    let mut saw_alt_ref = false;
    for variant in variants {
        match variant.gt {
            Some(0) => saw_ref_alt = true,
            Some(1) => saw_alt_ref = true,
            _ => {}
        }
    }

    if saw_ref_alt && saw_alt_ref {
        true
    } else {
        eprintln!("No informative numeric phased GT found: switching to unphased model");
        false
    }
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

    // Match scipy.optimize.minimize_scalar(method="bounded") used by AHO.
    // In sparse donors the MLE can be far below 0.001, so imposing a floor
    // changes the null likelihood and every downstream p-value.
    scipy_bounded_minimize(objective, 0.0, 1.0, 1e-5, 500)
}

fn depth_standardization(raw_depths: &[u32]) -> Result<(f64, f64)> {
    if raw_depths.is_empty() {
        anyhow::bail!("cannot standardize linear depth without observations");
    }

    // Welford's algorithm avoids cancellation when depths are large and tightly clustered.
    let mut center = 0.0;
    let mut squared_deviations = 0.0;
    for (index, &depth) in raw_depths.iter().enumerate() {
        let depth = depth as f64;
        let count = (index + 1) as f64;
        let delta = depth - center;
        center += delta / count;
        squared_deviations += delta * (depth - center);
    }

    if !center.is_finite() {
        anyhow::bail!("linear depth center is non-finite");
    }
    let variance = squared_deviations / raw_depths.len() as f64;
    let scale = variance.sqrt();
    if !scale.is_finite() {
        anyhow::bail!("linear depth scale is non-finite");
    }
    if scale <= 0.0 {
        anyhow::bail!("linear dispersion requires non-zero raw depth variance");
    }
    Ok((center, scale))
}

/// Optimize standardized linear dispersion parameters using Nelder-Mead.
///
/// Fits rho = expit(d1 + d2 * ((N - center) / scale)) using raw read depth for
/// the predictor and pseudocounted observations for the beta-binomial likelihood.
fn optimize_dispersion_linear(
    ref_counts: &[u32],
    n_array: &[u32],
    raw_depths: &[u32],
) -> Result<(f64, f64, f64, f64)> {
    if ref_counts.len() != n_array.len() || ref_counts.len() != raw_depths.len() {
        anyhow::bail!("linear dispersion inputs have inconsistent lengths");
    }
    let (depth_center, depth_scale) = depth_standardization(raw_depths)?;

    let objective = |params: [f64; 2]| -> f64 {
        let (d1, d2) = (params[0], params[1]);
        if !d1.is_finite() || !d2.is_finite() {
            return f64::INFINITY;
        }
        let mut neg_ll = 0.0;

        for ((&ref_count, &n), &raw_depth) in
            ref_counts.iter().zip(n_array.iter()).zip(raw_depths.iter())
        {
            let standardized_depth = (raw_depth as f64 - depth_center) / depth_scale;
            let predictor = d1 + standardized_depth * d2;
            if predictor.is_nan() {
                return f64::INFINITY;
            }
            let rho = clamp_rho(expit(predictor));
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

    // A constant-rho fit provides a data-driven intercept. Standardization makes
    // fixed slope starts meaningful across donors with different depth ranges.
    let scalar_rho = clamp_rho(optimize_dispersion(ref_counts, n_array)?);
    let scalar_intercept = (scalar_rho / (1.0 - scalar_rho)).ln();
    let starts: [[f64; 2]; 9] = [
        [scalar_intercept, 0.0],
        [scalar_intercept, -0.25],
        [scalar_intercept, 0.25],
        [scalar_intercept, -1.0],
        [scalar_intercept, 1.0],
        [scalar_intercept - 2.0, 0.0],
        [scalar_intercept + 2.0, 0.0],
        [0.0, -0.5],
        [0.0, 0.5],
    ];
    let mut best: Option<(f64, (f64, f64))> = None;
    for s in starts {
        if let Ok((d1, d2)) = nelder_mead_2d(objective, s, 1e-6, 1000) {
            let val = objective([d1, d2]);
            if val.is_finite() && best.is_none_or(|(bv, _)| val < bv) {
                best = Some((val, (d1, d2)));
            }
        }
    }
    let (d1, d2) = best
        .map(|(_, parameters)| parameters)
        .ok_or_else(|| anyhow::anyhow!("linear dispersion optimization failed from all starts"))?;
    Ok((d1, d2, depth_center, depth_scale))
}

/// Optimize probability parameter for alternative model
///
/// Python equivalent: `parse_opt()` calling `minimize_scalar(opt_prob, ...)`
fn optimize_prob(ref_counts: &[u32], n_array: &[u32], disp: f64) -> Result<(f64, f64)> {
    // For single SNP, optimize directly
    if ref_counts.len() == 1 {
        let objective = |prob: f64| -> f64 {
            opt_prob(prob, disp, ref_counts[0], n_array[0]).unwrap_or(f64::INFINITY)
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
/// * `disp` - Per-SNP dispersion (rho) slice, aligned 1:1 with ref/n/gt. A
///   length-1 slice broadcasts to all SNPs (matches the single/per-region scalar
///   path); a multi-element slice supplies a distinct rho per SNP (linear model's
///   per-row rho, mirroring Python's per-row `df["disp"]`).
///
/// # Returns
/// Tuple of (alt_ll, mu) - alternative-model log-likelihood and imbalance proportion
fn optimize_prob_phased(
    ref_counts: &[u32],
    n_array: &[u32],
    gt_array: &[u8],
    disp: &[f64],
) -> Result<(f64, f64)> {
    // Normalize phase so the first SNP is with respect to ref (matches Python).
    let flip = gt_array.first().is_some_and(|&g| g > 0);
    let gt_norm: Vec<u8> = gt_array
        .iter()
        .map(|&g| if flip { 1 - g } else { g })
        .collect();

    // objective(prob) = Σ opt_prob(|prob - gt_i|, rho_i, ref_i, n_i)
    // rho_i is the per-SNP dispersion (disp[i]) when a per-row slice is supplied,
    // or disp[0] broadcast across all SNPs for the scalar/single case.
    let objective = |prob: f64| -> f64 {
        let mut sum = 0.0;
        for (i, ((&k, &n), &g)) in ref_counts
            .iter()
            .zip(n_array.iter())
            .zip(gt_norm.iter())
            .enumerate()
        {
            let rho = if disp.len() > 1 { disp[i] } else { disp[0] };
            let phased_prob = (prob - g as f64).abs();
            match opt_prob(phased_prob, rho, k, n) {
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
/// This marginalizes over unknown phase independently at every position. Both
/// orientations receive equal prior weight, including the first SNP, so the
/// feature likelihood does not depend on input row order.
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
    if disp.is_empty() {
        anyhow::bail!("at least one dispersion value is required");
    }

    // Split per-variant dispersion for first vs subsequent positions (linear model)
    let (first_disp, phase_disp): (f64, &[f64]) = if disp.len() > 1 {
        (disp[0], &disp[1..])
    } else {
        // Scalar dispersion: use same value for all positions
        (disp[0], disp)
    };

    // Accumulate every position under both orientation hypotheses in log-space.
    // For each SNP:
    //   log(0.5 * p1 + 0.5 * p2) = ln(0.5) + logaddexp(ln(p1), ln(p2))
    let orientation_log_likelihood = |ref_count: u32, n: u32, rho: f64| -> Result<f64> {
        let rho = clamp_rho(rho);
        let alpha1 = prob * (1.0 - rho) / rho;
        let beta1 = (1.0 - prob) * (1.0 - rho) / rho;
        let bb1 = BetaBinomial::new(n, alpha1, beta1)
            .context("Failed to create beta-binomial for orientation 1")?;
        let lp1 = bb1.ln_f(&(ref_count as u64));

        let bb2 = BetaBinomial::new(n, beta1, alpha1)
            .context("Failed to create beta-binomial for orientation 2")?;
        let lp2 = bb2.ln_f(&(ref_count as u64));
        Ok((0.5_f64).ln() + logaddexp(lp1, lp2))
    };

    let mut log_likelihood = orientation_log_likelihood(first_ref, first_n, first_disp)?;

    for (i, (&ref_count, &n)) in phase_ref.iter().zip(phase_n.iter()).enumerate() {
        let rho = if phase_disp.len() > 1 {
            phase_disp[i.min(phase_disp.len() - 1)]
        } else {
            phase_disp[0]
        };
        log_likelihood += orientation_log_likelihood(ref_count, n, rho)?;
    }

    if !log_likelihood.is_finite() {
        return Ok(f64::INFINITY);
    }
    Ok(-log_likelihood)
}

/// Numerically stable two-argument log-sum-exp: ln(exp(a) + exp(b)).
///
/// Equivalent to NumPy's `np.logaddexp`. Subtracting the max before exponentiating
/// avoids overflow/underflow. If the max is non-finite (e.g. -inf when both
/// inputs are -inf), the max is returned directly.
#[inline]
fn logaddexp(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if !m.is_finite() {
        m
    } else {
        m + ((a - m).exp() + (b - m).exp()).ln()
    }
}

/// Exact two-sided beta-binomial p-value under the balanced null.
///
/// Outcomes at least as far from `n / 2` as the observed reference count are
/// included. Summation and normalization are performed in log-space so fitted
/// rho values close to the binomial boundary remain usable.
pub fn exact_beta_binomial_pvalue(k: u32, n: u32, rho: f64) -> Result<f64> {
    if k > n {
        anyhow::bail!("reference count cannot exceed total count");
    }

    let rho = clamp_rho(rho);
    let alpha = 0.5 * (1.0 - rho) / rho;
    let distribution = BetaBinomial::new(n, alpha, alpha)
        .context("Failed to create beta-binomial for exact SNV test")?;
    let observed_deviation = (2_i64 * k as i64 - n as i64).abs();
    let mut log_total = f64::NEG_INFINITY;
    let mut log_tail = f64::NEG_INFINITY;

    for outcome in 0..=n {
        let log_pmf = distribution.ln_f(&(outcome as u64));
        log_total = logaddexp(log_total, log_pmf);
        let deviation = (2_i64 * outcome as i64 - n as i64).abs();
        if deviation >= observed_deviation {
            log_tail = logaddexp(log_tail, log_pmf);
        }
    }

    if !log_total.is_finite() || !log_tail.is_finite() {
        anyhow::bail!("exact beta-binomial p-value was non-finite");
    }
    Ok((log_tail - log_total).exp().clamp(0.0, 1.0))
}

/// Bounded Brent minimizer matching SciPy's default scalar-bounded routine.
fn scipy_bounded_minimize<F>(f: F, x1: f64, x2: f64, xatol: f64, maxfun: usize) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    if !x1.is_finite() || !x2.is_finite() || x1 > x2 {
        anyhow::bail!("invalid bounded optimization interval");
    }
    let sqrt_eps = (2.2e-16_f64).sqrt();
    let golden_mean = 0.5 * (3.0 - 5.0_f64.sqrt());
    let (mut a, mut b) = (x1, x2);
    let mut fulc = a + golden_mean * (b - a);
    let (mut nfc, mut xf) = (fulc, fulc);
    let (mut rat, mut e) = (0.0_f64, 0.0_f64);
    let mut fx = f(xf);
    let (mut ffulc, mut fnfc) = (fx, fx);
    let mut fu = f64::INFINITY;
    let mut calls = 1_usize;
    let mut xm = 0.5 * (a + b);
    let mut tol1 = sqrt_eps * xf.abs() + xatol / 3.0;
    let mut tol2 = 2.0 * tol1;

    while (xf - xm).abs() > tol2 - 0.5 * (b - a) {
        let mut golden = true;
        if e.abs() > tol1 {
            let r = (xf - nfc) * (fx - ffulc);
            let q0 = (xf - fulc) * (fx - fnfc);
            let mut p = (xf - fulc) * q0 - (xf - nfc) * r;
            let mut q = 2.0 * (q0 - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let previous_e = e;
            e = rat;
            if p.abs() < (0.5 * q * previous_e).abs() && p > q * (a - xf) && p < q * (b - xf) {
                rat = p / q;
                let candidate = xf + rat;
                if candidate - a < tol2 || b - candidate < tol2 {
                    rat = tol1 * if xm - xf >= 0.0 { 1.0 } else { -1.0 };
                }
                golden = false;
            }
        }
        if golden {
            e = if xf >= xm { a - xf } else { b - xf };
            rat = golden_mean * e;
        }

        let direction = if rat >= 0.0 { 1.0 } else { -1.0 };
        let x = xf + direction * rat.abs().max(tol1);
        fu = f(x);
        calls += 1;
        if fu <= fx {
            if x >= xf {
                a = xf;
            } else {
                b = xf;
            }
            fulc = nfc;
            ffulc = fnfc;
            nfc = xf;
            fnfc = fx;
            xf = x;
            fx = fu;
        } else {
            if x < xf {
                a = x;
            } else {
                b = x;
            }
            if fu <= fnfc || nfc == xf {
                fulc = nfc;
                ffulc = fnfc;
                nfc = x;
                fnfc = fu;
            } else if fu <= ffulc || fulc == xf || fulc == nfc {
                fulc = x;
                ffulc = fu;
            }
        }

        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * xf.abs() + xatol / 3.0;
        tol2 = 2.0 * tol1;
        if calls >= maxfun {
            anyhow::bail!("bounded scalar optimization exceeded {maxfun} evaluations");
        }
    }
    if !xf.is_finite() || !fx.is_finite() || fu.is_nan() {
        anyhow::bail!("bounded scalar optimization produced a non-finite result");
    }
    Ok(xf)
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
            h *= inv_phi;
            c = a + inv_phi2 * h;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            h *= inv_phi;
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
        for point in simplex.iter_mut().skip(1) {
            let (x_i, _) = *point;
            let x_new = [
                x_best[0] + sigma * (x_i[0] - x_best[0]),
                x_best[1] + sigma * (x_i[1] - x_best[1]),
            ];
            *point = (x_new, f(x_new));
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
fn single_model_with_dispersion(
    variants: Vec<VariantCounts>,
    phased: bool,
) -> Result<(Vec<ImbalanceResult>, f64)> {
    if variants.is_empty() {
        return Ok((vec![], f64::NAN));
    }

    eprintln!("Optimizing dispersion parameter...");
    let parameters = fit_dispersion_parameters(&variants, AnalysisMethod::Single, 0)?;
    let DispersionParameters::Single { rho: disp } = parameters else {
        unreachable!("single fit returned linear parameters")
    };
    eprintln!("  Dispersion: {:.6}", disp);
    let results = test_fixed_dispersion(variants, phased, parameters, false, 0)?;
    Ok((results, disp))
}

pub fn single_model(variants: Vec<VariantCounts>, phased: bool) -> Result<Vec<ImbalanceResult>> {
    Ok(single_model_with_dispersion(variants, phased)?.0)
}

/// Linear dispersion model analysis
///
/// Python equivalent: `linear_model()` in as_analysis.py
/// Uses standardized raw depth for the globally optimized linear predictor.
pub fn linear_model_with_dispersion(
    variants: Vec<VariantCounts>,
    phased: bool,
) -> Result<(Vec<ImbalanceResult>, f64, f64)> {
    if variants.is_empty() {
        return Ok((vec![], f64::NAN, f64::NAN));
    }

    eprintln!("Optimizing linear dispersion parameters...");
    let parameters = fit_dispersion_parameters(&variants, AnalysisMethod::Linear, 0)?;
    let DispersionParameters::Linear { d1, d2, .. } = parameters else {
        unreachable!("linear fit returned scalar parameters")
    };
    eprintln!("  d1 = {:.6}, d2 = {:.6}", d1, d2);
    let results = test_fixed_dispersion(variants, phased, parameters, false, 0)?;
    Ok((results, d1, d2))
}

pub fn linear_model(variants: Vec<VariantCounts>, phased: bool) -> Result<Vec<ImbalanceResult>> {
    Ok(linear_model_with_dispersion(variants, phased)?.0)
}

fn fit_dispersion_parameters(
    variants: &[VariantCounts],
    method: AnalysisMethod,
    pseudocount: u32,
) -> Result<DispersionParameters> {
    if variants.is_empty() {
        anyhow::bail!("cannot fit dispersion without observations");
    }

    let ref_counts: Vec<u32> = variants.iter().map(|variant| variant.ref_count).collect();
    let n_array: Vec<u32> = variants
        .iter()
        .map(|variant| variant.ref_count + variant.alt_count)
        .collect();
    match method {
        AnalysisMethod::Single => Ok(DispersionParameters::Single {
            rho: optimize_dispersion(&ref_counts, &n_array)?,
        }),
        AnalysisMethod::Linear => {
            let raw_depths: Vec<u32> = variants
                .iter()
                .map(|variant| raw_depth(variant, pseudocount))
                .collect::<Result<_>>()?;
            let (d1, d2, depth_center, depth_scale) =
                optimize_dispersion_linear(&ref_counts, &n_array, &raw_depths)?;
            Ok(DispersionParameters::Linear {
                d1,
                d2,
                depth_center,
                depth_scale,
            })
        }
    }
}

fn raw_depth(variant: &VariantCounts, pseudocount: u32) -> Result<u32> {
    let raw_ref = variant.ref_count.checked_sub(pseudocount).ok_or_else(|| {
        anyhow::anyhow!("reference count is smaller than the configured pseudocount")
    })?;
    let raw_alt = variant.alt_count.checked_sub(pseudocount).ok_or_else(|| {
        anyhow::anyhow!("alternate count is smaller than the configured pseudocount")
    })?;
    raw_ref
        .checked_add(raw_alt)
        .ok_or_else(|| anyhow::anyhow!("raw read depth exceeds the supported range"))
}

fn canonicalize_region_indices(
    region_map: &mut HashMap<String, Vec<usize>>,
    variants: &[VariantCounts],
) {
    for indices in region_map.values_mut() {
        indices.sort_unstable_by(|left, right| {
            let left = &variants[*left];
            let right = &variants[*right];
            left.chrom
                .cmp(&right.chrom)
                .then_with(|| left.pos.cmp(&right.pos))
                .then_with(|| left.ref_count.cmp(&right.ref_count))
                .then_with(|| left.alt_count.cmp(&right.alt_count))
                .then_with(|| left.gt.cmp(&right.gt))
        });
    }
}

fn test_fixed_dispersion(
    variants: Vec<VariantCounts>,
    phased: bool,
    parameters: DispersionParameters,
    exact_snv_pvalues: bool,
    pseudocount: u32,
) -> Result<Vec<ImbalanceResult>> {
    if variants.is_empty() {
        return Ok(vec![]);
    }

    let ref_counts: Vec<u32> = variants.iter().map(|variant| variant.ref_count).collect();
    let n_array: Vec<u32> = variants
        .iter()
        .map(|variant| variant.ref_count + variant.alt_count)
        .collect();
    let raw_depths: Vec<u32> = variants
        .iter()
        .map(|variant| raw_depth(variant, pseudocount))
        .collect::<Result<_>>()?;
    let rho_array: Vec<f64> = raw_depths
        .iter()
        .map(|&depth| parameters.rho_at_depth(depth))
        .collect();

    let mut region_map: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, variant) in variants.iter().enumerate() {
        region_map
            .entry(variant.region.clone())
            .or_default()
            .push(i);
    }
    canonicalize_region_indices(&mut region_map, &variants);

    eprintln!(
        "Optimizing imbalance likelihood for {} regions...",
        region_map.len()
    );

    let results: Result<Vec<_>> = region_map
        .par_iter()
        .map(|(region, indices)| -> Result<ImbalanceResult> {
            let region_ref: Vec<u32> = indices.iter().map(|&i| ref_counts[i]).collect();
            let region_n: Vec<u32> = indices.iter().map(|&i| n_array[i]).collect();
            let region_rho: Vec<f64> = indices.iter().map(|&i| rho_array[i]).collect();

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

            let all_gt = indices.iter().all(|&i| variants[i].gt.is_some());
            let (alt_ll, mu) = if phased && indices.len() > 1 && all_gt {
                let region_gt: Vec<u8> = indices
                    .iter()
                    .map(|&i| variants[i].gt.expect("gt present (checked by all_gt)"))
                    .collect();
                optimize_prob_phased(&region_ref, &region_n, &region_gt, &region_rho)?
            } else if indices.len() > 1 {
                let first_ref = region_ref[0];
                let first_n = region_n[0];
                let phase_ref = &region_ref[1..];
                let phase_n = &region_n[1..];
                let objective = |prob: f64| -> f64 {
                    optimize_prob_unphased_dp(
                        prob,
                        &region_rho,
                        first_ref,
                        first_n,
                        phase_ref,
                        phase_n,
                    )
                    .unwrap_or(f64::INFINITY)
                };
                // The orientation-marginal likelihood is symmetric around 0.5.
                // Restrict to one half to make the reported magnitude identifiable.
                let mu = golden_section_search(objective, 0.5, 1.0, 1e-6)?;
                let alt_ll = -objective(mu);
                (alt_ll, mu)
            } else {
                optimize_prob(&region_ref, &region_n, region_rho[0])?
            };

            let lrt = (-2.0 * (null_ll - alt_ll)).max(0.0);
            let pval = if exact_snv_pvalues && indices.len() == 1 {
                let variant = &variants[indices[0]];
                let raw_ref = variant.ref_count.checked_sub(pseudocount).ok_or_else(|| {
                    anyhow::anyhow!("reference count is smaller than the pseudocount")
                })?;
                let raw_n = raw_depth(variant, pseudocount)?;
                exact_beta_binomial_pvalue(raw_ref, raw_n, parameters.rho_at_depth(raw_n))?
            } else {
                let chi2 =
                    ChiSquared::new(1.0).context("Failed to create chi-squared distribution")?;
                1.0 - chi2.cdf(lrt)
            };

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
    results.sort_by(|left, right| left.region.cmp(&right.region));
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
struct PreparedAnalysis {
    variants: Vec<VariantCounts>,
    effective_phased: bool,
}

fn prepare_analysis(variants: Vec<VariantCounts>, config: &AnalysisConfig) -> PreparedAnalysis {
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

    let phased = effective_phased(&filtered, config.phased);

    // Deduplicate variants for Python parity.
    // Python's get_imbalance() does: df[keep_cols].drop_duplicates()
    let before = filtered.len();
    let deduped = deduplicate_variants(filtered, phased);
    if deduped.len() < before {
        eprintln!(
            "Deduplicated {} → {} variants ({} duplicates removed)",
            before,
            deduped.len(),
            before - deduped.len()
        );
    }
    PreparedAnalysis {
        variants: deduped,
        effective_phased: phased,
    }
}

fn remove_result_pseudocounts(results: &mut [ImbalanceResult], pseudocount: u32) {
    for result in results {
        let region_pseudocount = pseudocount.saturating_mul(result.snp_count as u32);
        if result.ref_count < region_pseudocount
            || result.alt_count < region_pseudocount
            || result.n < 2 * region_pseudocount
        {
            eprintln!(
                "[WARN] Counts smaller than pseudocount for region {}: ref={}, alt={}, n={}, pc={}",
                result.region, result.ref_count, result.alt_count, result.n, region_pseudocount
            );
        }
        result.ref_count = result.ref_count.saturating_sub(region_pseudocount);
        result.alt_count = result.alt_count.saturating_sub(region_pseudocount);
        result.n = result.n.saturating_sub(2 * region_pseudocount);
    }
}

fn analyze_prepared(
    prepared: PreparedAnalysis,
    config: &AnalysisConfig,
    parameters: DispersionParameters,
    exact_snv_pvalues: bool,
) -> Result<AnalysisOutput> {
    let n_observations = prepared.variants.len();
    let mut results = test_fixed_dispersion(
        prepared.variants,
        prepared.effective_phased,
        parameters,
        exact_snv_pvalues,
        config.pseudocount,
    )?;
    remove_result_pseudocounts(&mut results, config.pseudocount);

    Ok(AnalysisOutput {
        results,
        rho: parameters.rho(),
        linear_params: parameters.linear_params(),
        linear_depth_center: parameters
            .linear_depth_standardization()
            .map(|standardization| standardization.0),
        linear_depth_scale: parameters
            .linear_depth_standardization()
            .map(|standardization| standardization.1),
        n_observations,
        effective_phased: prepared.effective_phased,
    })
}

/// Fit only the balanced-null nuisance parameters for one independent table.
pub fn fit_imbalance_dispersion(
    variants: Vec<VariantCounts>,
    config: &AnalysisConfig,
) -> Result<DispersionFit> {
    let prepared = prepare_analysis(variants, config);
    if prepared.variants.is_empty() {
        anyhow::bail!("no variants passed the count threshold for dispersion fitting");
    }

    match config.method {
        AnalysisMethod::Single => eprintln!("Optimizing dispersion parameter..."),
        AnalysisMethod::Linear => eprintln!("Optimizing linear dispersion parameters..."),
    }
    let parameters =
        fit_dispersion_parameters(&prepared.variants, config.method, config.pseudocount)?;
    match parameters {
        DispersionParameters::Single { rho } => eprintln!("  Dispersion: {:.6}", rho),
        DispersionParameters::Linear {
            d1,
            d2,
            depth_center,
            depth_scale,
        } => {
            eprintln!(
                "  d1 = {:.6}, d2 = {:.6}, depth center = {:.6}, depth scale = {:.6}",
                d1, d2, depth_center, depth_scale
            )
        }
    }
    Ok(DispersionFit {
        parameters,
        n_observations: prepared.variants.len(),
    })
}

/// Test effects with nuisance parameters supplied by a separate fitting run.
pub fn analyze_imbalance_with_fixed_dispersion(
    variants: Vec<VariantCounts>,
    config: &AnalysisConfig,
    parameters: DispersionParameters,
    exact_snv_pvalues: bool,
) -> Result<AnalysisOutput> {
    let parameters = parameters.validate_for(config.method)?;
    let prepared = prepare_analysis(variants, config);
    analyze_prepared(prepared, config, parameters, exact_snv_pvalues)
}

pub fn analyze_imbalance_detailed(
    variants: Vec<VariantCounts>,
    config: &AnalysisConfig,
) -> Result<AnalysisOutput> {
    let prepared = prepare_analysis(variants, config);
    let n_observations = prepared.variants.len();
    if n_observations == 0 {
        return Ok(AnalysisOutput {
            results: vec![],
            rho: None,
            linear_params: None,
            linear_depth_center: None,
            linear_depth_scale: None,
            n_observations,
            effective_phased: prepared.effective_phased,
        });
    }

    match config.method {
        AnalysisMethod::Single => eprintln!("Optimizing dispersion parameter..."),
        AnalysisMethod::Linear => eprintln!("Optimizing linear dispersion parameters..."),
    }
    let parameters =
        fit_dispersion_parameters(&prepared.variants, config.method, config.pseudocount)?;
    match parameters {
        DispersionParameters::Single { rho } => eprintln!("  Dispersion: {:.6}", rho),
        DispersionParameters::Linear {
            d1,
            d2,
            depth_center,
            depth_scale,
        } => {
            eprintln!(
                "  d1 = {:.6}, d2 = {:.6}, depth center = {:.6}, depth scale = {:.6}",
                d1, d2, depth_center, depth_scale
            )
        }
    }
    analyze_prepared(prepared, config, parameters, false)
}

pub fn analyze_imbalance(
    variants: Vec<VariantCounts>,
    config: &AnalysisConfig,
) -> Result<Vec<ImbalanceResult>> {
    Ok(analyze_imbalance_detailed(variants, config)?.results)
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
    fn test_standardized_linear_predictor_is_algebraically_equivalent() {
        let raw_intercept = -4.0;
        let raw_slope = 0.02;
        let depth_center = 20.0;
        let depth_scale = 5.0;
        let parameters = DispersionParameters::Linear {
            d1: raw_intercept + raw_slope * depth_center,
            d2: raw_slope * depth_scale,
            depth_center,
            depth_scale,
        };

        for depth in [0, 10, 20, 40, 100] {
            let expected = clamp_rho(expit(raw_intercept + raw_slope * depth as f64));
            assert!((parameters.rho_at_depth(depth) - expected).abs() < 1e-15);
        }
    }

    #[test]
    fn test_linear_predictor_is_not_clipped() {
        let parameters = DispersionParameters::Linear {
            d1: 0.0,
            d2: 100.0,
            depth_center: 0.0,
            depth_scale: 1.0,
        };
        assert!(parameters.rho_at_depth(1) > 0.999999);
    }

    #[test]
    fn test_depth_standardization_rejects_zero_variance() {
        let error = depth_standardization(&[20, 20, 20]).unwrap_err();
        assert!(error.to_string().contains("non-zero raw depth variance"));
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
    fn test_scipy_bounded_minimize_preserves_near_zero_solution() {
        let minimum = scipy_bounded_minimize(|rho| rho, 0.0, 1.0, 1e-5, 500).unwrap();
        let scipy_boundary = 5.9608609865491405e-6;
        assert!((minimum - scipy_boundary).abs() < 1e-15);
        assert!(minimum < 0.001);
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

    #[test]
    fn test_exact_beta_binomial_pvalue_is_symmetric_and_handles_low_rho() {
        let left = exact_beta_binomial_pvalue(2, 10, 0.08).unwrap();
        let right = exact_beta_binomial_pvalue(8, 10, 0.08).unwrap();
        assert!((left - right).abs() < 1e-12);
        assert!((exact_beta_binomial_pvalue(5, 10, 0.08).unwrap() - 1.0).abs() < 1e-12);

        // At the clamp boundary this converges to the ordinary two-sided
        // binomial probability for observing 0 or 10 successes.
        let near_binomial = exact_beta_binomial_pvalue(0, 10, 0.0).unwrap();
        assert!((near_binomial - 2.0 / 1024.0).abs() < 1e-7);
    }

    fn vc_peak(pos: u32, rc: u32, ac: u32, gt: Option<u8>) -> VariantCounts {
        VariantCounts {
            chrom: "chr1".to_string(),
            pos,
            ref_count: rc,
            alt_count: ac,
            region: "peak1".to_string(),
            gt,
        }
    }

    fn assert_result_close(a: &ImbalanceResult, b: &ImbalanceResult) {
        assert_eq!(a.region, b.region);
        assert_eq!(a.ref_count, b.ref_count);
        assert_eq!(a.alt_count, b.alt_count);
        assert_eq!(a.n, b.n);
        assert!((a.null_ll - b.null_ll).abs() < 1e-10);
        assert!((a.alt_ll - b.alt_ll).abs() < 1e-10);
        assert!((a.mu - b.mu).abs() < 1e-10);
        assert!((a.lrt - b.lrt).abs() < 1e-10);
        assert!((a.pval - b.pval).abs() < 1e-10);
    }

    #[test]
    fn test_effective_phased_requires_both_numeric_orientations() {
        let both = vec![vc_peak(100, 18, 4, Some(0)), vc_peak(200, 5, 15, Some(1))];
        assert!(effective_phased(&both, true));
        assert!(!effective_phased(&both, false));

        let none = vec![vc_peak(100, 18, 4, None), vc_peak(200, 5, 15, None)];
        assert!(!effective_phased(&none, true));

        let same = vec![vc_peak(100, 18, 4, Some(0)), vc_peak(200, 12, 8, Some(0))];
        assert!(!effective_phased(&same, true));

        let partial = vec![vc_peak(100, 18, 4, Some(0)), vc_peak(200, 12, 8, None)];
        assert!(!effective_phased(&partial, true));
    }

    #[test]
    fn test_analyze_imbalance_all_same_numeric_gt_falls_back_to_unphased() {
        let variants = vec![vc_peak(100, 18, 4, Some(0)), vc_peak(200, 12, 8, Some(0))];
        let phased_cfg = AnalysisConfig {
            min_count: 0,
            pseudocount: 0,
            method: AnalysisMethod::Single,
            phased: true,
        };
        let unphased_cfg = AnalysisConfig {
            phased: false,
            ..phased_cfg.clone()
        };

        let phased = analyze_imbalance(variants.clone(), &phased_cfg).unwrap();
        let unphased = analyze_imbalance(variants, &unphased_cfg).unwrap();
        assert_eq!(phased.len(), 1);
        assert_eq!(unphased.len(), 1);
        assert_result_close(&phased[0], &unphased[0]);
    }

    #[test]
    fn test_single_model_missing_gt_uses_unphased_multi_snp_path() {
        let variants = vec![vc_peak(100, 18, 4, None), vc_peak(200, 5, 15, None)];
        let phased = single_model(variants.clone(), true).unwrap();
        let unphased = single_model(variants, false).unwrap();
        assert_eq!(phased.len(), 1);
        assert_eq!(unphased.len(), 1);
        assert_result_close(&phased[0], &unphased[0]);
    }

    #[test]
    fn test_multi_snp_output_removes_every_pseudocount() {
        let variants = vec![vc_peak(100, 10, 5, None), vc_peak(200, 7, 8, None)];
        let config = AnalysisConfig {
            min_count: 0,
            pseudocount: 1,
            method: AnalysisMethod::Single,
            phased: false,
        };
        let output = analyze_imbalance_detailed(variants, &config).unwrap();
        assert_eq!(output.results.len(), 1);
        let result = &output.results[0];
        assert_eq!(result.snp_count, 2);
        assert_eq!(result.ref_count, 17);
        assert_eq!(result.alt_count, 13);
        assert_eq!(result.n, 30);
    }

    #[test]
    fn test_unphased_feature_is_permutation_invariant_with_linear_nuisance() {
        let variants = vec![
            vc_peak(300, 17, 3, None),
            vc_peak(100, 4, 16, None),
            vc_peak(200, 13, 12, None),
        ];
        let mut reversed = variants.clone();
        reversed.reverse();
        let config = AnalysisConfig {
            min_count: 0,
            pseudocount: 0,
            method: AnalysisMethod::Linear,
            phased: false,
        };
        let nuisance = DispersionParameters::Linear {
            d1: -3.0,
            d2: 0.025,
            depth_center: 20.0,
            depth_scale: 5.0,
        };

        let forward =
            analyze_imbalance_with_fixed_dispersion(variants, &config, nuisance, false).unwrap();
        let backward =
            analyze_imbalance_with_fixed_dispersion(reversed, &config, nuisance, false).unwrap();
        assert_eq!(forward.results.len(), 1);
        assert_eq!(backward.results.len(), 1);
        assert_result_close(&forward.results[0], &backward.results[0]);
        assert_eq!(forward.results[0].fdr_pval, backward.results[0].fdr_pval);
    }

    #[test]
    fn test_fixed_nuisance_exact_snv_uses_raw_counts() {
        let variants = vec![vc_peak(100, 8, 2, None)];
        let no_pseudocount = AnalysisConfig {
            min_count: 0,
            pseudocount: 0,
            method: AnalysisMethod::Single,
            phased: false,
        };
        let with_pseudocount = AnalysisConfig {
            pseudocount: 3,
            ..no_pseudocount.clone()
        };
        let nuisance = DispersionParameters::Single { rho: 0.04 };

        let raw = analyze_imbalance_with_fixed_dispersion(
            variants.clone(),
            &no_pseudocount,
            nuisance,
            true,
        )
        .unwrap();
        let padded =
            analyze_imbalance_with_fixed_dispersion(variants, &with_pseudocount, nuisance, true)
                .unwrap();
        assert_eq!(raw.results[0].ref_count, 8);
        assert_eq!(padded.results[0].ref_count, 8);
        assert_eq!(raw.results[0].alt_count, 2);
        assert_eq!(padded.results[0].alt_count, 2);
        assert!((raw.results[0].pval - padded.results[0].pval).abs() < 1e-12);
        assert_eq!(raw.rho, Some(0.04));
        assert_eq!(padded.rho, Some(0.04));
    }

    #[test]
    fn test_fixed_linear_nuisance_is_reported_without_refitting() {
        let config = AnalysisConfig {
            min_count: 0,
            pseudocount: 0,
            method: AnalysisMethod::Linear,
            phased: false,
        };
        let nuisance = DispersionParameters::Linear {
            d1: -4.25,
            d2: 0.031,
            depth_center: 10.0,
            depth_scale: 2.0,
        };
        let output = analyze_imbalance_with_fixed_dispersion(
            vec![vc_peak(100, 7, 5, None)],
            &config,
            nuisance,
            false,
        )
        .unwrap();
        assert_eq!(output.linear_params, Some((-4.25, 0.031)));
        assert_eq!(output.linear_depth_center, Some(10.0));
        assert_eq!(output.linear_depth_scale, Some(2.0));
        assert_eq!(output.rho, None);
        assert_eq!(output.n_observations, 1);
    }

    #[test]
    fn test_fixed_linear_nuisance_rejects_invalid_scale() {
        let config = AnalysisConfig {
            min_count: 0,
            pseudocount: 0,
            method: AnalysisMethod::Linear,
            phased: false,
        };
        for depth_scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let error = analyze_imbalance_with_fixed_dispersion(
                vec![vc_peak(100, 7, 5, None)],
                &config,
                DispersionParameters::Linear {
                    d1: -4.0,
                    d2: 0.1,
                    depth_center: 12.0,
                    depth_scale,
                },
                false,
            )
            .unwrap_err();
            assert!(error.to_string().contains("finite and positive"));
        }
    }

    #[test]
    fn test_linear_fit_reports_raw_depth_standardization() {
        let variants = vec![
            vc_peak(100, 8, 2, None),
            vc_peak(200, 10, 10, None),
            vc_peak(300, 18, 12, None),
        ];
        let config = AnalysisConfig {
            min_count: 0,
            pseudocount: 3,
            method: AnalysisMethod::Linear,
            phased: false,
        };
        let fit = fit_imbalance_dispersion(variants, &config).unwrap();
        let (center, scale) = fit.parameters.linear_depth_standardization().unwrap();
        assert!((center - 20.0).abs() < 1e-12);
        assert!((scale - (200.0_f64 / 3.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_linear_fit_rejects_zero_raw_depth_variance() {
        let variants = vec![vc_peak(100, 8, 2, None), vc_peak(200, 6, 4, None)];
        let config = AnalysisConfig {
            min_count: 0,
            pseudocount: 2,
            method: AnalysisMethod::Linear,
            phased: false,
        };
        let error = fit_imbalance_dispersion(variants, &config).unwrap_err();
        assert!(error.to_string().contains("non-zero raw depth variance"));
    }
}
