//! Positive-domain beta-binomial arithmetic from SciPy 1.13.1.
//!
//! Sources: scipy/special/special/cephes/{gamma,beta,polevl}.h and
//! scipy/stats/_discrete_distns.py at tag v1.13.1. Preserve coefficients,
//! branch boundaries, and arithmetic order: algebraic rearrangement can
//! change which of two symmetric optima a bounded optimizer selects.
//!
//! Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
//! Cephes Gamma: Copyright 1984, 1987, 1989, 1992 Stephen L. Moshier.
//! Cephes Beta: Copyright 1984, 1987 Stephen L. Moshier.
//! See the distributed third-party SciPy license for redistribution terms.

// These literals intentionally preserve the upstream decimal coefficients.
#![allow(clippy::excessive_precision)]

use anyhow::{bail, Result};

const GAMMA_P: [f64; 7] = [
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1,
];
const GAMMA_Q: [f64; 8] = [
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0,
];
const GAMMA_STIR: [f64; 5] = [
    7.87311395793093628397E-4,
    -2.29549961613378126380E-4,
    -2.68132617805781232825E-3,
    3.47222221605458667310E-3,
    8.33333333333482257126E-2,
];
const GAMMA_A: [f64; 5] = [
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2,
];
const GAMMA_B: [f64; 6] = [
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5,
];
const GAMMA_C: [f64; 6] = [
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6,
];
const MAXGAM: f64 = 171.624376956302725;
// Retain the pinned Cephes literal on Rust versions without EULER_GAMMA.
#[allow(clippy::approx_constant)]
const CEPHES_EULER_GAMMA: f64 = 0.5772156649015329;

/// Reproduce the measured reference compiler's contraction policy. The pinned
/// manylinux x86_64 SciPy wheel evaluates Cephes operations separately, whereas
/// the macOS arm64 reference contracts these expressions. This is independent
/// of the dataset and preserves the original arithmetic on other targets.
#[inline]
fn cephes_mul_add(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    {
        a * b + c
    }
    #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
    {
        a.mul_add(b, c)
    }
}

fn polevl(x: f64, coefficients: &[f64]) -> f64 {
    let mut result = coefficients[0];
    for &coefficient in &coefficients[1..] {
        result = cephes_mul_add(result, x, coefficient);
    }
    result
}

fn p1evl(x: f64, coefficients: &[f64]) -> f64 {
    let mut result = x + coefficients[0];
    for &coefficient in &coefficients[1..] {
        result = cephes_mul_add(result, x, coefficient);
    }
    result
}

fn stirf(x: f64) -> f64 {
    if x >= MAXGAM {
        return f64::INFINITY;
    }
    let mut w = 1.0 / x;
    w = cephes_mul_add(w, polevl(w, &GAMMA_STIR), 1.0);
    let mut y = x.exp();
    if x > 143.01608 {
        let v = x.powf(0.5 * x - 0.25);
        y = v * (v / y);
    } else {
        y = x.powf(x - 0.5) / y;
    }
    2.50662827463100050242 * y * w
}

fn gamma_positive(mut x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    if x > 33.0 {
        return stirf(x);
    }
    let mut z = 1.0;
    while x >= 3.0 {
        x -= 1.0;
        z *= x;
    }
    while x < 2.0 {
        if x < 1e-9 {
            return z / (cephes_mul_add(CEPHES_EULER_GAMMA, x, 1.0) * x);
        }
        z /= x;
        x += 1.0;
    }
    if x == 2.0 {
        return z;
    }
    x -= 2.0;
    z * polevl(x, &GAMMA_P) / polevl(x, &GAMMA_Q)
}

fn lgamma_positive(mut x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    if x < 13.0 {
        let mut z = 1.0;
        let mut p = 0.0;
        let mut u = x;
        while u >= 3.0 {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while u < 2.0 {
            if u == 0.0 {
                return f64::INFINITY;
            }
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if u == 2.0 {
            return z.ln();
        }
        p -= 2.0;
        x += p;
        p = x * polevl(x, &GAMMA_B) / p1evl(x, &GAMMA_C);
        return z.ln() + p;
    }
    if x > 2.556348e305 {
        return f64::INFINITY;
    }
    let mut q = cephes_mul_add(x - 0.5, x.ln(), -x) + 0.91893853320467274178;
    if x > 1.0e8 {
        return q;
    }
    let p = 1.0 / (x * x);
    if x >= 1000.0 {
        q += ((7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3) * p
            + 0.0833333333333333333333)
            / x;
    } else {
        q += polevl(p, &GAMMA_A) / x;
    }
    q
}

// Preserve the upstream operand order rather than rewriting as +=.
#[allow(clippy::assign_op_pattern)]
fn betaln_positive(mut a: f64, mut b: f64) -> f64 {
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }
    if a > 1e6 * b && a > 1e6 {
        let mut r = lgamma_positive(b);
        r -= b * a.ln();
        r += b * (1.0 - b) / (2.0 * a);
        r += b * (1.0 - b) * (1.0 - 2.0 * b) / (12.0 * a * a);
        r += -b * b * (1.0 - b) * (1.0 - b) / (12.0 * a * a * a);
        return r;
    }
    let mut y = a + b;
    if y > MAXGAM || a > MAXGAM || b > MAXGAM {
        y = lgamma_positive(y);
        y = lgamma_positive(b) - y;
        y = lgamma_positive(a) + y;
        return y;
    }
    y = gamma_positive(y);
    a = gamma_positive(a);
    b = gamma_positive(b);
    if y == 0.0 {
        return f64::INFINITY;
    }
    if (a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs() {
        y = b / y;
        y *= a;
    } else {
        y = a / y;
        y *= b;
    }
    y.abs().ln()
}

/// NumPy 1.26.4 contiguous float64 pairwise reduction (128-value blocks).
/// Source: numpy/core/src/umath/loops_utils.h.src; see NUMPY-LICENSE.txt.
pub(crate) fn numpy_sum(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 8 {
        return values.iter().fold(-0.0, |sum, &value| sum + value);
    }
    if n <= 128 {
        let mut lanes = [0.0; 8];
        lanes.copy_from_slice(&values[..8]);
        let end = n - n % 8;
        for block in values[8..end].as_chunks::<8>().0 {
            for (lane, &value) in lanes.iter_mut().zip(block) {
                *lane += value;
            }
        }
        let mut result = ((lanes[0] + lanes[1]) + (lanes[2] + lanes[3]))
            + ((lanes[4] + lanes[5]) + (lanes[6] + lanes[7]));
        for &value in &values[end..] {
            result += value;
        }
        return result;
    }
    let mut middle = n / 2;
    middle -= middle % 8;
    numpy_sum(&values[..middle]) + numpy_sum(&values[middle..])
}

/// SciPy's beta-binomial log PMF, preserving its three betaln evaluations.
pub(crate) fn betabinom_logpmf(k: u32, n: u32, alpha: f64, beta: f64) -> Result<f64> {
    if !alpha.is_finite() || !beta.is_finite() || alpha <= 0.0 || beta <= 0.0 {
        bail!("beta-binomial shape parameters must be finite and positive");
    }
    if k > n {
        return Ok(f64::NEG_INFINITY);
    }
    if n == 0 {
        return Ok(0.0);
    }
    let k = f64::from(k);
    let n = f64::from(n);
    let combiln = -(n + 1.0).ln() - betaln_positive(n - k + 1.0, k + 1.0);
    Ok(combiln + betaln_positive(k + alpha, n - k + beta) - betaln_positive(alpha, beta))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 5e-14 * expected.abs().max(1.0));
    }

    #[test]
    fn scipy_1131_gamma_and_loggamma_branches() {
        // Generated with scipy.special.gamma/gammaln in SciPy 1.13.1.
        for (x, gamma, lgamma) in [
            (1e-12, 999999999999.4226, 27.63102111592797),
            (0.1, 9.513507698668732, 2.252712651734206),
            (2.5, 1.329340388179137, 0.2846828704729192),
            (13.0, 479001600.0, 19.98721449566189),
            (33.01, 2.724528906472356e35, 81.59277378960284),
            (143.02, 2.976421070788908e245, 565.2240693793439),
            (171.0, 7.257415615308e306, 706.5730622457875),
        ] {
            close(gamma_positive(x), gamma);
            close(lgamma_positive(x), lgamma);
        }
        close(lgamma_positive(1000.0), 5905.220423209181);
        close(lgamma_positive(1e8), 1742068066.1038349);
    }

    #[test]
    fn scipy_1131_beta_branches() {
        for (a, b, expected) in [
            (0.1, 0.2, 2.680978479346915),
            (30.0, 40.0, -48.30174909591612),
            (80.0, 100.0, -124.629853209889),
            (1e8, 0.2, -2.160072325559689),
            (1e-8, 10.0, 18.42068071566268),
            (500.0, 500.0, -694.988722485713),
        ] {
            close(betaln_positive(a, b), expected);
            close(betaln_positive(b, a), expected);
        }
    }

    #[test]
    fn numpy_1264_reduction_boundaries() {
        // Independently generated np.sum(float64) references, including cancellation.
        for (n, expected) in [
            (7, 1e16),
            (8, 0.0),
            (16, 1e16),
            (127, 1.0000000000000004e16),
            (128, 0.0),
            (129, 0.0),
            (256, 1e16),
            (502, 1.0000000000000004e16),
            (513, 0.0),
        ] {
            let values: Vec<_> = (0..n)
                .map(|i| match i % 3 {
                    0 => 1e16,
                    1 => -1e16,
                    _ => 1.0,
                })
                .collect();
            assert_eq!(numpy_sum(&values), expected);
        }
    }

    #[test]
    fn real_probe_components() {
        let a = 8.977092566062371;
        let b = 5.5481483259805575;
        let x = 8.0 + a;
        let y = 4.0 + b;
        eprintln!(
            "betaln={:?}, {:?}, {:?}; gamma={:?}, {:?}, {:?}; numerator shapes={:?},{:?}",
            betaln_positive(5.0, 9.0),
            betaln_positive(x, y),
            betaln_positive(a, b),
            gamma_positive(x),
            gamma_positive(y),
            gamma_positive(x + y),
            x,
            y
        );
        // Independently measured SciPy 1.13.1 values in both reference builds.
        let expected = if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
            (19621650968466.156, 132624.78469324368)
        } else {
            (19621650968466.152, 132624.78469324365)
        };
        assert_eq!(gamma_positive(x), expected.0);
        assert_eq!(gamma_positive(y), expected.1);
        close(betaln_positive(x, y), -17.307163194505247);
    }

    #[test]
    fn scipy_1131_beta_binomial_expression() {
        for (k, n, a, b, expected) in [
            (0, 10, 0.1, 0.2, -0.656035869515514),
            (5, 10, 1.0, 1.0, -2.3978952727983707),
            (100, 200, 50.0, 50.0, -3.4271728015926755),
            (12, 50, 1e7, 1e7, -9.134989546611905),
            (3, 12, 0.0001, 0.5, -10.168874588765625),
        ] {
            close(betabinom_logpmf(k, n, a, b).unwrap(), expected);
        }
    }

    #[test]
    fn beta_binomial_support_and_parameter_validation() {
        assert_eq!(betabinom_logpmf(0, 0, 1.0, 2.0).unwrap(), 0.0);
        assert_eq!(betabinom_logpmf(2, 1, 1.0, 2.0).unwrap(), f64::NEG_INFINITY);
        for shape in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(betabinom_logpmf(0, 1, shape, 1.0).is_err());
            assert!(betabinom_logpmf(0, 1, 1.0, shape).is_err());
        }
    }
}
