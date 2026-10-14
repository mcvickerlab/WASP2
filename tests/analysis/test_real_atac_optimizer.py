"""Conditional bulk inference against archived Python on real ATAC observations.

The independently frozen Python fixture fixes the nuisance dispersion to isolate
allelic-proportion optimization. Full-data comparisons separately retain strict
mu checks: this test does not certify orientation for symmetric unphased peaks.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scripts.compare_archive_analysis import compare_frames

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
    if case["family"] == "balanced_anchor" and not case["phased"]:
        # A balanced anchor makes the unphased objective symmetric. Check all
        # statistics and either equally likely optimum, without choosing a side.
        retained = counts.loc[counts.ref_count + counts.alt_count >= 10]
        first = retained.groupby("region", sort=False).first()
        assert (first.ref_count == first.alt_count).all()
        comparison = compare_frames(expected.drop(columns="mu"), actual.drop(columns="mu"))
        expected_mu = expected.set_index("region").mu
        actual_mu = actual.set_index("region").mu.reindex(expected_mu.index)
        same = np.isclose(actual_mu, expected_mu, rtol=5e-5, atol=5e-5)
        mirrored = np.isclose(actual_mu, 1 - expected_mu, rtol=5e-5, atol=5e-5)
        assert (same | mirrored).all()
    else:
        comparison = compare_frames(expected, actual)
    assert comparison["status"] == "PASS", comparison
