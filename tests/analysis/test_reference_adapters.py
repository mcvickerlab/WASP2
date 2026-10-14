"""Shape repairs must preserve the original probability recurrence exactly."""

import numpy as np
import pytest
from scripts.compare_archive_analysis import adapt_linear_reference, load_archive


@pytest.mark.parametrize("mu", [0.06, 0.25, 0.38196601125, 0.7, 0.94])
@pytest.mark.parametrize("rho", [0.03, 0.99])
def test_scalar_shape_adapter_preserves_archive_arithmetic(mu, rho):
    modules = load_archive()
    original = modules["as_analysis"].opt_unphased_dp
    first_ref, first_n = np.array([10]), np.array([20])
    rest_ref, rest_n = np.array([3, 18, 0]), np.array([20, 24, 10])
    expected = original(mu, rho, first_ref, first_n, rest_ref, rest_n)
    adapt_linear_reference(modules)
    actual = modules["as_analysis"].opt_unphased_dp(mu, rho, first_ref, first_n, rest_ref, rest_n)
    assert actual == expected


@pytest.mark.parametrize("mu", [0.06, 0.25, 0.38196601125, 0.7, 0.94])
def test_vector_adapter_only_slices_dispersion(mu):
    modules = load_archive()
    probability = modules["as_analysis"].opt_prob
    dispersions = np.array([0.01, 0.03, 0.1, 0.99])
    first_ref, first_n = np.array([10]), np.array([20])
    rest_ref, rest_n = np.array([3, 18, 0]), np.array([20, 24, 10])
    first = probability(mu, dispersions[0], first_ref[0], first_n[0])
    forward = probability(mu, dispersions[1:], rest_ref, rest_n, log=False)
    reverse = probability(1 - mu, dispersions[1:], rest_ref, rest_n, log=False)
    previous = 1
    for p1, p2 in zip(forward, reverse, strict=True):
        phase1, phase2 = previous * p1, previous * p2
        previous = 0.5 * phase1 + 0.5 * phase2
    expected = float(first - np.log(previous))
    adapt_linear_reference(modules)
    actual = modules["as_analysis"].opt_unphased_dp(
        mu, dispersions, first_ref, first_n, rest_ref, rest_n
    )
    assert actual == expected
