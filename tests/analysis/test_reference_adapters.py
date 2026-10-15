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


@pytest.mark.parametrize("mu", [0.06, 0.25, 0.38196601125, 0.7, 0.94])
def test_python_unphased_likelihood_preserves_archive_recurrence(mu):
    from analysis.as_analysis import opt_unphased_dp

    archive = load_archive()["as_analysis"]
    ref, depth = np.array([10, 3, 18, 0]), np.array([20, 20, 24, 10])
    expected = archive.opt_unphased_dp(mu, 0.03, ref[:1], depth[:1], ref[1:], depth[1:])
    actual = opt_unphased_dp(mu, 0.03, ref[:1], depth[:1], ref[1:], depth[1:])
    assert actual == float(np.asarray(expected).item())


def test_python_unphased_underflow_fallback_remains_finite():
    from scipy.stats import betabinom

    from analysis.as_analysis import opt_unphased_dp

    mu, rho = 0.2, 0.03
    ref, depth = np.arange(2000, 2200), np.full(200, 100_000)
    shape = (1 - rho) / rho
    forward = betabinom.logpmf(ref, depth, mu * shape, (1 - mu) * shape)
    reverse = betabinom.logpmf(ref, depth, (1 - mu) * shape, mu * shape)
    expected = -forward[0] - np.sum(np.logaddexp(forward[1:], reverse[1:]) - np.log(2))
    actual = opt_unphased_dp(mu, rho, ref[:1], depth[:1], ref[1:], depth[1:])
    assert np.isfinite(actual)
    assert np.isclose(actual, expected, rtol=5e-5, atol=5e-5)


def test_wide_archive_accumulators_prevent_overflow(monkeypatch):
    from scipy.sparse import csr_matrix
    from scripts.compare_archive_analysis import single_cell_counts

    archive = load_archive(widen_single_cell_sums=True)["as_analysis_sc"]
    historical = single_cell_counts()
    historical.layers["ref"] = csr_matrix(np.full(historical.shape, 40000, dtype=np.uint16))

    def inspect_fit(function, args, **kwargs):
        ref, depth = args
        np.testing.assert_array_equal(ref, np.full(historical.n_obs, 160001))
        assert ref.dtype == np.uint64
        assert np.all(depth >= ref)
        raise RuntimeError("inspected reference inputs")

    monkeypatch.setattr(archive, "minimize_scalar", inspect_fit)
    with pytest.raises(RuntimeError, match="inspected reference inputs"):
        archive.get_imbalance_sc(historical, sample="sample1", groups=["A", "B"], pseudocount=1)
