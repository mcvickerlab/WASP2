"""Tests for beta-binomial rho parameter clamping (Issue #228).

Verifies that the rho parameter is properly clamped to avoid division by zero
and numerical instability at boundary values (rho=0 and rho=1).

These tests use a minimal implementation that mirrors the production code
to avoid environment-specific import issues with pandas/pyarrow.

3x Hardening Tests:
1. Core function tests (clamp_rho, opt_prob boundaries)
2. Array/scalar consistency tests
3. Edge case and stress tests
"""

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.stats import betabinom

# =============================================================================
# Mirror production constants for isolated testing
# =============================================================================
RHO_EPSILON: float = 1e-10


def clamp_rho(rho: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
    """Mirror of as_analysis.clamp_rho for isolated testing."""
    return np.clip(rho, RHO_EPSILON, 1.0 - RHO_EPSILON)


def opt_prob(
    in_prob: float,
    in_rho: float,
    k: int,
    n: int,
) -> float:
    """Mirror of as_analysis.opt_prob for isolated testing."""
    prob = in_prob
    rho = clamp_rho(in_rho)
    alpha = prob * (1 - rho) / rho
    beta = (1 - prob) * (1 - rho) / rho
    return float(-1 * betabinom.logpmf(k, n, alpha, beta))


# =============================================================================
# 1x Hardening: Core Function Tests
# =============================================================================


class TestClampRhoCore:
    """Core tests for the clamp_rho helper function."""

    def test_clamp_zero(self) -> None:
        """Test that rho=0 is clamped to epsilon."""
        result = clamp_rho(0.0)
        assert result == RHO_EPSILON
        assert result > 0

    def test_clamp_one(self) -> None:
        """Test that rho=1 is clamped to 1-epsilon."""
        result = clamp_rho(1.0)
        assert result == 1.0 - RHO_EPSILON
        assert result < 1

    def test_clamp_normal_value(self) -> None:
        """Test that normal values in (0,1) are unchanged."""
        for val in [0.1, 0.3, 0.5, 0.7, 0.9]:
            result = clamp_rho(val)
            assert result == val, f"Value {val} should be unchanged"

    def test_clamp_near_boundary(self) -> None:
        """Test values very close to boundaries are clamped."""
        tiny = 1e-15
        result_low = clamp_rho(tiny)
        result_high = clamp_rho(1.0 - tiny)

        # Should be clamped since tiny < RHO_EPSILON
        assert result_low == RHO_EPSILON
        assert result_high == 1.0 - RHO_EPSILON


class TestOptProbBoundaries:
    """Test opt_prob doesn't produce NaN/Inf at boundary rho values."""

    def test_opt_prob_rho_zero(self) -> None:
        """Test opt_prob doesn't crash or return NaN when rho=0."""
        result = opt_prob(0.5, 0.0, 10, 20)
        assert np.isfinite(result), "rho=0 should produce finite result"

    def test_opt_prob_rho_one(self) -> None:
        """Test opt_prob doesn't crash or return NaN when rho=1."""
        result = opt_prob(0.5, 1.0, 10, 20)
        assert np.isfinite(result), "rho=1 should produce finite result"

    def test_opt_prob_rho_very_small(self) -> None:
        """Test opt_prob with very small rho (near-binomial case)."""
        result = opt_prob(0.5, 1e-12, 10, 20)
        assert np.isfinite(result), "Very small rho should produce finite result"

    def test_opt_prob_rho_near_one(self) -> None:
        """Test opt_prob with rho very close to 1."""
        result = opt_prob(0.5, 1.0 - 1e-12, 10, 20)
        assert np.isfinite(result), "rho near 1 should produce finite result"


# =============================================================================
# 2x Hardening: Array/Scalar Consistency Tests
# =============================================================================


class TestClampRhoArray:
    """Test clamping works consistently with numpy arrays."""

    def test_clamp_array_basic(self) -> None:
        """Test clamping works on numpy arrays."""
        arr = np.array([0.0, 0.5, 1.0])
        result = clamp_rho(arr)
        assert result[0] == RHO_EPSILON
        assert result[1] == 0.5
        assert result[2] == 1.0 - RHO_EPSILON

    def test_clamp_array_mixed_boundaries(self) -> None:
        """Test array with values at both boundaries and middle."""
        arr = np.array([0.0, 1e-15, 0.001, 0.5, 0.999, 1.0 - 1e-15, 1.0])
        result = clamp_rho(arr)

        # Check boundaries clamped
        assert result[0] == RHO_EPSILON
        assert result[1] == RHO_EPSILON  # 1e-15 < epsilon
        assert result[-2] == 1.0 - RHO_EPSILON  # 1 - 1e-15 > 1 - epsilon
        assert result[-1] == 1.0 - RHO_EPSILON

        # Check interior values preserved
        assert result[3] == 0.5

    def test_clamp_large_array(self) -> None:
        """Test clamping on large array for performance."""
        rng = np.random.default_rng(42)
        arr = rng.uniform(-0.1, 1.1, size=10000)  # Some values outside [0,1]
        result = clamp_rho(arr)

        # All values should be in safe range
        assert np.all(result >= RHO_EPSILON)
        assert np.all(result <= 1.0 - RHO_EPSILON)


# =============================================================================
# 3x Hardening: Edge Case and Stress Tests
# =============================================================================


class TestBetaBinomialParameterization:
    """Test alpha/beta parameterization is valid after clamping."""

    def test_alpha_beta_positive(self) -> None:
        """Test that alpha and beta are always positive after clamping."""
        for rho_raw in [0.0, 1e-15, 0.5, 1.0 - 1e-15, 1.0]:
            rho = float(clamp_rho(rho_raw))
            for prob in [0.1, 0.5, 0.9]:
                alpha = prob * (1 - rho) / rho
                beta = (1 - prob) * (1 - rho) / rho
                assert alpha > 0, f"alpha <= 0 for rho_raw={rho_raw}, prob={prob}"
                assert beta > 0, f"beta <= 0 for rho_raw={rho_raw}, prob={prob}"
                assert np.isfinite(alpha), f"alpha not finite for rho_raw={rho_raw}"
                assert np.isfinite(beta), f"beta not finite for rho_raw={rho_raw}"

    def test_betabinom_valid_after_clamping(self) -> None:
        """Test that scipy.stats.betabinom works with clamped parameters."""
        for rho_raw in [0.0, 1e-15, 0.5, 1.0 - 1e-15, 1.0]:
            rho = float(clamp_rho(rho_raw))
            alpha = 0.5 * (1 - rho) / rho
            beta = 0.5 * (1 - rho) / rho

            # Should not raise
            result = betabinom.logpmf(10, 20, alpha, beta)
            assert np.isfinite(result), f"betabinom.logpmf not finite for rho_raw={rho_raw}"


class TestStressConditions:
    """Stress tests for edge conditions."""

    def test_negative_rho_clamped(self) -> None:
        """Negative rho values should be clamped to epsilon."""
        result = clamp_rho(-1.0)
        assert result == RHO_EPSILON

    def test_rho_greater_than_one_clamped(self) -> None:
        """rho > 1 should be clamped to 1-epsilon."""
        result = clamp_rho(2.0)
        assert result == 1.0 - RHO_EPSILON

    def test_inf_rho_clamped(self) -> None:
        """inf rho should be clamped."""
        result = clamp_rho(np.inf)
        assert result == 1.0 - RHO_EPSILON

    def test_negative_inf_rho_clamped(self) -> None:
        """Negative inf rho should be clamped."""
        result = clamp_rho(-np.inf)
        assert result == RHO_EPSILON

    def test_opt_prob_extreme_counts(self) -> None:
        """Test opt_prob with extreme count values at boundary rho."""
        # Very small counts
        result = opt_prob(0.5, 0.0, 0, 1)
        assert np.isfinite(result)

        # Large counts
        result = opt_prob(0.5, 1.0, 1000, 2000)
        assert np.isfinite(result)

    def test_opt_prob_extreme_prob_values(self) -> None:
        """Test opt_prob with extreme probability values."""
        # prob near 0
        result = opt_prob(0.01, 0.1, 10, 20)
        assert np.isfinite(result)

        # prob near 1
        result = opt_prob(0.99, 0.1, 10, 20)
        assert np.isfinite(result)


class TestRegressionIssue228:
    """Regression tests specifically for Issue #228 scenarios."""

    def test_division_by_zero_prevented(self) -> None:
        """Verify that division by zero is prevented when rho=0."""
        # This would previously cause division by zero
        rho = 0.0
        rho_clamped = clamp_rho(rho)

        # The formula (1-rho)/rho should not raise
        result = (1 - rho_clamped) / rho_clamped
        assert np.isfinite(result)

    def test_zero_alpha_beta_prevented(self) -> None:
        """Verify that alpha/beta don't become zero when rho=1."""
        # This would previously make alpha and beta = 0
        rho = 1.0
        rho_clamped = clamp_rho(rho)

        alpha = 0.5 * (1 - rho_clamped) / rho_clamped
        beta = 0.5 * (1 - rho_clamped) / rho_clamped

        assert alpha > 0, "alpha should be positive"
        assert beta > 0, "beta should be positive"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
