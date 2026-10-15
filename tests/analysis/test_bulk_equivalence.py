"""Verification must distinguish symmetric optima from arbitrary mu mismatches."""

import numpy as np
import pandas as pd
import pytest
from scripts.compare_archive_analysis import (
    compare_bulk_equivalence,
    prepare_bulk_reference_counts,
)


def example():
    expected = pd.DataFrame(
        [
            {
                "region": "peak",
                "ref_count": 12,
                "alt_count": 20,
                "N": 32,
                "snp_count": 2,
                "mu": 0.3,
                "null_ll": -5.0,
                "alt_ll": -4.0,
                "lrt": 2.0,
                "pval": 0.157,
                "fdr_pval": 0.157,
            }
        ]
    )
    observed = expected.copy()
    observed["mu"] = 0.7
    counts = pd.DataFrame(
        [
            {"chrom": "chr21", "pos": 1, "region": "peak", "ref_count": 7, "alt_count": 7},
            {"chrom": "chr21", "pos": 2, "region": "peak", "ref_count": 5, "alt_count": 13},
        ]
    )
    retained = prepare_bulk_reference_counts(
        counts, min_count=10, pseudocount=0, effective_phased=False, region_col="region"
    )
    return expected, observed, retained


def verify(expected, observed, retained, **kwargs):
    return compare_bulk_equivalence(
        expected,
        observed,
        retained_counts=retained,
        effective_phased=kwargs.pop("effective_phased", False),
        reference_log_likelihood=kwargs.pop("objective", lambda region, mu: -4.0),
        dispersion_comparison=kwargs.pop("dispersion", {"status": "PASS"}),
        **kwargs,
    )


def test_mirrored_optimum_keeps_strict_failure_and_records_evidence():
    expected, observed, retained = example()
    before = observed.copy(deep=True)
    result = verify(expected, observed, retained)
    assert result["strict_comparison"]["status"] == "FAIL"
    assert result["strict_comparison"]["failures"] == ["mu"]
    assert result["statistical_equivalence"] == "PASS"
    evidence = result["regions"][0]
    assert evidence["status"] == "SYMMETRIC_EQUIVALENT"
    assert evidence["reference_log_likelihoods"] == [-4.0, -4.0, -4.0]
    pd.testing.assert_frame_equal(observed, before)


@pytest.mark.parametrize(
    "failure",
    [
        "asymmetric",
        "phased",
        "not_complementary",
        "nonfinite_mu",
        "wrong_likelihood",
        "nonfinite_objective",
        "incorrect_reported_likelihood",
        "counts",
        "dispersion",
        "missing_anchor",
        "missing_mu",
        "different_regions",
        "wrong_observation_count",
        "invalid_retained_depth",
    ],
)
def test_invalid_equivalence_is_rejected(failure):
    expected, observed, retained = example()
    kwargs = {}
    if failure == "asymmetric":
        retained.loc[0, "ref_count"] += 1
    elif failure == "phased":
        kwargs["effective_phased"] = True
    elif failure == "not_complementary":
        observed["mu"] = 0.8
    elif failure == "nonfinite_mu":
        observed["mu"] = np.nan
    elif failure == "wrong_likelihood":
        kwargs["objective"] = lambda region, mu: -4.0 + mu
    elif failure == "nonfinite_objective":
        kwargs["objective"] = lambda region, mu: np.nan
    elif failure == "incorrect_reported_likelihood":
        kwargs["objective"] = lambda region, mu: -8.0
    elif failure == "counts":
        observed["ref_count"] += 1
    elif failure == "dispersion":
        kwargs["dispersion"] = {"status": "FAIL"}
    elif failure == "missing_anchor":
        retained = retained.drop(columns="ref_count")
    elif failure == "missing_mu":
        observed = observed.drop(columns="mu")
    elif failure == "different_regions":
        observed["region"] = "other"
    elif failure == "wrong_observation_count":
        retained = pd.concat([retained, retained.iloc[[1]]])
    elif failure == "invalid_retained_depth":
        retained.loc[0, "N"] += 1
    assert verify(expected, observed, retained, **kwargs)["statistical_equivalence"] == "FAIL"


def test_matching_mu_stays_strict_without_symmetry_evidence():
    expected, _, retained = example()
    retained.loc[0, "ref_count"] += 1
    result = verify(expected, expected.copy(), retained, effective_phased=True)
    assert result["strict_comparison"]["status"] == "PASS"
    assert result["statistical_equivalence"] == "PASS"
    assert result["regions"][0]["status"] == "STRICT_MATCH"


@pytest.mark.parametrize("pc", [0, 1])
@pytest.mark.parametrize("phased", [False, True])
def test_retained_anchor_is_selected_after_filtering_and_phase_aware_deduplication(pc, phased):
    counts = pd.DataFrame(
        [
            {
                "chrom": "chr21",
                "pos": 10,
                "region": "peak",
                "ref_count": 1,
                "alt_count": 2,
                "GT": "0|1",
            },
            {
                "chrom": "chr21",
                "pos": 20,
                "region": "peak",
                "ref_count": 5,
                "alt_count": 5,
                "GT": "0|1",
            },
            {
                "chrom": "chr21",
                "pos": 20,
                "region": "peak",
                "ref_count": 5,
                "alt_count": 5,
                "GT": "0|1",
            },
            {
                "chrom": "chr21",
                "pos": 20,
                "region": "peak",
                "ref_count": 5,
                "alt_count": 5,
                "GT": "1|0",
            },
            {
                "chrom": "chr21",
                "pos": 30,
                "region": "peak",
                "ref_count": 3,
                "alt_count": 8,
                "GT": "0|1",
            },
        ]
    )
    retained = prepare_bulk_reference_counts(
        counts, min_count=10, pseudocount=pc, effective_phased=phased, region_col="region"
    )
    assert retained.iloc[0].pos == 20
    assert retained.iloc[0].ref_count == retained.iloc[0].alt_count == 5 + pc
    assert retained.pos.tolist() == ([20, 20, 30] if phased else [20, 30])
    assert "GT" in retained if phased else "GT" not in retained
