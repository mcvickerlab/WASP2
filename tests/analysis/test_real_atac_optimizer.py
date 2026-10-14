"""Conditional bulk inference against archived Python on real ATAC observations.

The independently frozen Python fixture fixes the nuisance dispersion to isolate
allelic-proportion optimization. Full-data comparisons separately retain strict
mu checks: this test does not certify orientation for symmetric unphased peaks.
"""

import json
from pathlib import Path

import pandas as pd
import pytest
from scipy.special import expit
from scripts.compare_archive_analysis import (
    adapt_linear_reference,
    compare_bulk_equivalence,
    load_archive,
    prepare_bulk_reference_counts,
)

REFERENCE = json.loads(
    (Path(__file__).parent / "fixtures/real_atac_optimizer_reference.json").read_text()
)


@pytest.mark.rust
@pytest.mark.parametrize("case", REFERENCE["cases"], ids=lambda case: case["name"])
def test_real_atac_conditional_optimizer_matches_python(tmp_path, case):
    import wasp2_rust

    counts = pd.DataFrame(case["counts"])
    count_file = tmp_path / "counts.tsv"
    counts.to_csv(count_file, sep="\t", index=False)
    output = wasp2_rust.analyze_imbalance_run(
        str(count_file),
        method=case["method"],
        phased=case["phased"],
        pseudocount=case["pseudocount"],
        region_col="region",
        **case["fixed_dispersion"],
    )
    actual = pd.DataFrame(output["results"])
    expected = pd.DataFrame(case["expected"])
    retained = prepare_bulk_reference_counts(
        counts,
        min_count=10,
        pseudocount=case["pseudocount"],
        effective_phased=output["effective_phased"],
        region_col="region",
    )
    modules = load_archive()
    if case["method"] == "linear":
        adapt_linear_reference(modules)
    archive = modules["as_analysis"]

    def reference_likelihood(region, mu):
        rows = retained.loc[retained.region == region]
        ref = rows.ref_count.to_numpy()
        depth = rows.N.to_numpy()
        fixed = case["fixed_dispersion"]
        rho = (
            fixed["rho"]
            if case["method"] == "single"
            else expit(
                fixed["linear_d1"]
                + fixed["linear_d2"]
                * (
                    (depth - 2 * case["pseudocount"] - fixed["linear_depth_center"])
                    / fixed["linear_depth_scale"]
                )
            )
        )
        if len(rows) == 1:
            value = archive.opt_prob(mu, rho, ref[0], depth[0])
        else:
            value = archive.opt_unphased_dp(mu, rho, ref[:1], depth[:1], ref[1:], depth[1:])
        return -float(value)

    # Fixed nuisance parameters are identical by construction in this conditional
    # test; the full-data runner separately verifies fitted dispersion predictions.
    comparison = compare_bulk_equivalence(
        expected,
        actual,
        retained_counts=retained,
        effective_phased=output["effective_phased"],
        reference_log_likelihood=reference_likelihood,
        dispersion_comparison={"status": "PASS", "source": "identical fixed nuisance"},
    )
    assert comparison["statistical_equivalence"] == "PASS", comparison
    if comparison["strict_comparison"]["status"] != "PASS":
        assert comparison["strict_comparison"]["failures"] == ["mu"]
        assert any(row["status"] == "SYMMETRIC_EQUIVALENT" for row in comparison["regions"])
