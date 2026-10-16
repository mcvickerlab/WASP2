"""Strict Python compatibility and independent numerical-stability checks.

Conditional expectations execute archived Python in the current environment.
The original frozen Mac expectations remain unchanged and are additionally
checked in their pinned environment; neither oracle is generated from Rust.
The high-depth checks use SciPy's log PMF as an independent mathematical oracle;
they do not claim compatibility with the archive's underflowing product loop.
"""

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy
from scipy.stats import betabinom
from scripts.compare_archive_analysis import compare_frames, fixed_bulk_reference

REFERENCE = json.loads(
    (Path(__file__).parent / "fixtures/real_atac_optimizer_reference.json").read_text()
)


@pytest.mark.rust
@pytest.mark.parametrize("case", REFERENCE["cases"], ids=lambda case: case["name"])
def test_fixed_real_atac_matches_native_python_and_frozen_mac(tmp_path, case):
    import wasp2_rust

    counts_file = tmp_path / "counts.tsv"
    pd.DataFrame(case["counts"]).to_csv(counts_file, sep="\t", index=False)
    output = wasp2_rust.analyze_imbalance_run(
        str(counts_file),
        method=case["method"],
        phased=case["phased"],
        pseudocount=case["pseudocount"],
        region_col="region",
        **case["fixed_dispersion"],
    )
    expected = fixed_bulk_reference(case)
    comparison = compare_frames(expected, pd.DataFrame(output["results"]))
    assert comparison["status"] == "PASS", comparison
    # Preserve the original frozen Mac expectations without applying its
    # compiler/library bit patterns to another platform's archived Python.
    if (platform.system(), platform.machine(), np.__version__, scipy.__version__) == (
        "Darwin",
        "arm64",
        "1.26.4",
        "1.13.1",
    ):
        frozen = compare_frames(pd.DataFrame(case["expected"]), pd.DataFrame(output["results"]))
        assert frozen["status"] == "PASS", frozen


@pytest.mark.rust
@pytest.mark.parametrize("rho", [1e-8, 0.1, 0.99])
@pytest.mark.parametrize("ref_count,total", [(0, 10), (10, 10), (1, 1000), (500, 1000)])
def test_fixed_dispersion_null_logpmf_matches_scipy(tmp_path, rho, ref_count, total):
    import wasp2_rust

    counts_file = tmp_path / "counts.tsv"
    pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "pos": 101,
                "ref": "A",
                "alt": "G",
                "region": "peak",
                "ref_count": ref_count,
                "alt_count": total - ref_count,
            }
        ]
    ).to_csv(counts_file, sep="\t", index=False)
    output = wasp2_rust.analyze_imbalance_run(
        str(counts_file),
        method="single",
        phased=False,
        pseudocount=0,
        region_col="region",
        rho=rho,
    )
    actual = output["results"][0]["null_ll"]
    shape = 0.5 * (1.0 - rho) / rho
    expected = betabinom.logpmf(ref_count, total, shape, shape)
    assert np.isfinite(actual)
    assert np.isclose(actual, expected, rtol=5e-5, atol=5e-5)


@pytest.mark.rust
def test_many_high_depth_unphased_observations_remain_finite(tmp_path):
    import wasp2_rust

    # Unique positions ensure these are distinct retained observations. A long
    # product underflows even though the equivalent log likelihood is finite.
    total = 100_000
    ref_counts = np.arange(200, dtype=np.int64) + 2000
    counts_file = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "chrom": "chr1",
            "pos": np.arange(200) + 101,
            "ref": "A",
            "alt": "G",
            "region": "peak",
            "ref_count": ref_counts,
            "alt_count": total - ref_counts,
        }
    ).to_csv(counts_file, sep="\t", index=False)
    rho = 0.1
    output = wasp2_rust.analyze_imbalance_run(
        str(counts_file),
        method="single",
        phased=False,
        pseudocount=0,
        region_col="region",
        rho=rho,
    )
    result = output["results"][0]
    assert result["snp_count"] == len(ref_counts)
    assert all(np.isfinite(result[key]) for key in ["mu", "null_ll", "alt_ll", "pval"])
    assert 0.0 < result["mu"] < 1.0
    shape = 0.5 * (1.0 - rho) / rho
    expected_null = np.sum(betabinom.logpmf(ref_counts, total, shape, shape))
    assert np.isclose(result["null_ll"], expected_null, rtol=5e-5, atol=5e-5)

    alpha = result["mu"] * (1.0 - rho) / rho
    beta = (1.0 - result["mu"]) * (1.0 - rho) / rho
    forward = betabinom.logpmf(ref_counts, total, alpha, beta)
    reverse = betabinom.logpmf(ref_counts, total, beta, alpha)
    expected_alt = forward[0] + np.sum(np.logaddexp(forward[1:], reverse[1:]) - np.log(2.0))
    assert np.isclose(result["alt_ll"], expected_alt, rtol=5e-5, atol=5e-5)
    product = 1.0
    for p1, p2 in zip(np.exp(forward[1:]), np.exp(reverse[1:]), strict=True):
        product = 0.5 * (product * p1) + 0.5 * (product * p2)
    assert product == 0.0


def decimal_beta_binomial_logpmf(k, n, alpha, beta):
    """Independent rising-factorial PMF; no SciPy or gamma/beta function reuse."""
    from decimal import Decimal, localcontext
    from math import comb

    with localcontext() as context:
        context.prec = 70
        a, b = Decimal(str(alpha)), Decimal(str(beta))
        value = Decimal(comb(n, k))
        for i in range(k):
            value *= a + i
        for i in range(n - k):
            value *= b + i
        for i in range(n):
            value /= a + b + i
        return float(value.ln())


@pytest.mark.rust
@pytest.mark.parametrize("ref_count,total", [(0, 0), (0, 14), (4, 8), (7, 14), (14, 14)])
def test_likelihood_precision_against_independent_decimal(tmp_path, ref_count, total):
    import wasp2_rust
    from analysis.as_analysis import opt_prob

    rho, alpha = 0.1, 4.5
    expected = decimal_beta_binomial_logpmf(ref_count, total, alpha, alpha)
    vector = opt_prob(
        0.5, rho, np.array([ref_count], dtype=np.uint64), np.array([total], dtype=np.uint64)
    )
    scalar = opt_prob(0.5, rho, np.uint64(ref_count), np.uint64(total))
    assert float(np.asarray(vector).item()) == scalar
    assert scalar == pytest.approx(-expected, abs=1e-12, rel=1e-12)
    path = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "chrom": ["chr1"],
            "pos": [101],
            "ref": ["A"],
            "alt": ["G"],
            "region": ["peak"],
            "ref_count": [ref_count],
            "alt_count": [total - ref_count],
        }
    ).to_csv(path, sep="\t", index=False)
    result = wasp2_rust.analyze_imbalance_run(
        str(path), method="single", min_count=0, pseudocount=0, region_col="region", rho=rho
    )["results"][0]
    assert result["null_ll"] == pytest.approx(expected, abs=1e-12, rel=1e-12)
    if 2 * ref_count == total:
        assert result["alt_ll"] == pytest.approx(result["null_ll"], abs=1e-12, rel=1e-12)
        assert result["pval"] == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("first_ref,symmetric", [(3, True), (4, False)])
def test_unphased_anchor_symmetry_is_conditional(first_ref, symmetric):
    from analysis.as_analysis import opt_unphased_dp

    rho = 0.1
    ref = np.array([first_ref, 1, 0, 0, 0], dtype=np.uint64)
    depth = np.array([6, 20, 0, 0, 1], dtype=np.uint64)
    values = []
    for mu in (0.2, 0.8):
        expected_first = decimal_beta_binomial_logpmf(
            int(ref[0]), int(depth[0]), mu * 9, (1 - mu) * 9
        )
        rest = [
            np.logaddexp(
                decimal_beta_binomial_logpmf(int(k), int(n), mu * 9, (1 - mu) * 9),
                decimal_beta_binomial_logpmf(int(k), int(n), (1 - mu) * 9, mu * 9),
            )
            - np.log(2)
            for k, n in zip(ref[1:], depth[1:], strict=True)
        ]
        expected = -expected_first - sum(rest)
        actual = opt_unphased_dp(mu, rho, ref[:1], depth[:1], ref[1:], depth[1:])
        assert actual == pytest.approx(expected, abs=1e-12, rel=1e-12)
        values.append(actual)
    if symmetric:
        assert values[0] == pytest.approx(values[1], abs=1e-12, rel=1e-12)
    else:
        assert abs(values[0] - values[1]) > 1.0


@pytest.mark.skipif(
    (np.__version__, scipy.__version__) != ("1.26.4", "1.13.1"),
    reason="Historical dtype dispatch reproduction uses pinned NumPy/SciPy",
)
def test_untouched_archive_balanced_precision_defect_is_reported(tmp_path):
    from scripts.compare_archive_analysis import load_archive

    path = tmp_path / "balanced.tsv"
    pd.DataFrame(
        {
            "chrom": ["chr1"],
            "pos": [101],
            "ref": ["A"],
            "alt": ["G"],
            "region": ["peak"],
            "ref_count": [7],
            "alt_count": [7],
            "other_count": [0],
        }
    ).to_csv(path, sep="\t", index=False)
    narrow = load_archive()["as_analysis"].get_imbalance(
        str(path), region_col="region", min_count=0, pseudocount=0
    )
    wide = load_archive(widen_bulk_csv_counts=True)["as_analysis"].get_imbalance(
        str(path), region_col="region", min_count=0, pseudocount=0
    )
    assert float(narrow.pval.iloc[0]) < 0.9999
    assert float(wide.pval.iloc[0]) == 1.0
    assert compare_frames(narrow, wide)["status"] == "FAIL"


@pytest.mark.rust
@pytest.mark.parametrize("ref_count,alt_count,pc", [(40000, 40000, 1), (65536, 65537, 0)])
def test_wide_bulk_counts_preserve_raw_outputs_and_fitting_depth(
    tmp_path, ref_count, alt_count, pc
):
    import wasp2_rust

    path = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "chrom": ["chr1"],
            "pos": [101],
            "ref": ["A"],
            "alt": ["G"],
            "region": ["peak"],
            "ref_count": [ref_count],
            "alt_count": [alt_count],
        }
    ).to_csv(path, sep="\t", index=False)
    rho = 0.1
    result = wasp2_rust.analyze_imbalance_run(
        str(path), method="single", min_count=0, pseudocount=pc, region_col="region", rho=rho
    )["results"][0]
    assert (result["ref_count"], result["alt_count"], result["N"]) == (
        ref_count,
        alt_count,
        ref_count + alt_count,
    )
    shape = 0.5 * (1 - rho) / rho
    expected_null = betabinom.logpmf(
        np.uint64(ref_count + pc), np.uint64(ref_count + alt_count + 2 * pc), shape, shape
    )
    assert result["null_ll"] == pytest.approx(expected_null, abs=1e-9, rel=1e-9)
