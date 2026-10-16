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


@pytest.mark.parametrize("column", ["num_snps_group1", "num_snps_group2"])
def test_differential_snp_membership_counts_require_exact_agreement(column):
    import pandas as pd
    from scripts.compare_archive_analysis import compare_frames

    expected = pd.DataFrame({"region": ["feature"], column: [2], "fdr_pval": [0.5]})
    observed = expected.copy()
    observed[column] = 3
    comparison = compare_frames(expected, observed)
    assert comparison["status"] == "FAIL"
    assert column in comparison["failures"]


@pytest.mark.parametrize("column", ["ref_count", "N", "num_snps_group1", "num_snps_group2"])
@pytest.mark.parametrize("value", [2**53, 2**63, 2**64 - 2])
def test_count_verifier_preserves_uint64_precision(column, value):
    import pandas as pd
    from scripts.compare_archive_analysis import compare_frames

    expected = pd.DataFrame(
        {"region": ["feature"], column: pd.Series([value], dtype="uint64"), "fdr_pval": [0.5]}
    )
    observed = expected.copy()
    assert compare_frames(expected, observed)["status"] == "PASS"
    observed[column] = pd.Series([value + 1], dtype="uint64")
    result = compare_frames(expected, observed)
    assert result["status"] == "FAIL"
    assert result["max_absolute_differences"][column] == 1


@pytest.mark.parametrize("wide", [False, True])
def test_bulk_csv_adapter_only_changes_count_reader_dtypes(tmp_path, monkeypatch, wide):
    import pandas as pd

    path = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "chrom": ["chr1"],
            "pos": [101],
            "ref": ["A"],
            "alt": ["G"],
            "ref_count": [40000],
            "alt_count": [40000],
            "other_count": [0],
        }
    ).to_csv(path, sep="\t", index=False)
    module = load_archive(widen_bulk_csv_counts=wide)["as_analysis"]

    def inspect_fit(frame, *args, **kwargs):
        expected_type = np.dtype("uint64" if wide else "uint16")
        assert frame.ref_count.dtype == frame.alt_count.dtype == frame.N.dtype == expected_type
        assert int(frame.ref_count.iloc[0]) == 40001
        assert int(frame.N.iloc[0]) == (80002 if wide else 14466)
        raise RuntimeError("inspected bulk CSV fitting inputs")

    monkeypatch.setattr(module, "single_model", inspect_fit)
    with pytest.raises(RuntimeError, match="inspected bulk CSV"):
        module.get_imbalance(str(path), pseudocount=1)
    # Other modules are untouched by this independently selected adapter.
    other = load_archive(widen_bulk_csv_counts=wide)["as_analysis_sc"]
    from scripts.compare_archive_analysis import single_cell_counts

    def inspect_sc(function, args, **kwargs):
        assert all(value.dtype == np.uint16 for value in args)
        raise RuntimeError("single-cell sums untouched")

    monkeypatch.setattr(other, "minimize_scalar", inspect_sc)
    with pytest.raises(RuntimeError, match="single-cell sums untouched"):
        other.get_imbalance_sc(single_cell_counts(), sample="sample1", groups=["A", "B"])


def test_reference_classification_keeps_strict_failures():
    import pandas as pd
    from scripts.compare_archive_analysis import classify_reference_differences

    legacy = pd.DataFrame({"region": ["feature"], "mu": [0.2], "pval": [0.8], "fdr_pval": [0.8]})
    wide = legacy.assign(mu=0.8)
    result = classify_reference_differences(legacy, wide, wide.copy())
    assert result["legacy_strict_comparison"]["status"] == "FAIL"
    assert result["widened_comparison"]["status"] == "PASS"
    assert result["unexplained_rows"] == 0
    assert result["mismatched_rows"][0]["complementary_estimate_candidates"] == ["mu"]
    damaged = wide.assign(pval=0.3)
    assert classify_reference_differences(legacy, wide, damaged)["unexplained_rows"] == 1
