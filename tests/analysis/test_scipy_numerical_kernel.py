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
