"""Donor inference must match separate archived Python invocations.

Use two identical donor count tables so a shared dispersion fit also has the
same optimum as separate fits. References come from Python, not from Rust.
"""

import json
from pathlib import Path

import pandas as pd
import pytest
from scripts.compare_archive_analysis import bulk_counts, compare_frames

from analysis.run_analysis_per_donor import run_per_donor_analysis

REFERENCE = json.loads(
    (Path(__file__).parent / "fixtures/archive_python_reference.json").read_text()
)


@pytest.mark.rust
@pytest.mark.parametrize("unit,phased", [("snv", False), ("feature", False), ("feature", True)])
@pytest.mark.parametrize("dispersion_scope", ["global", "per-donor"])
@pytest.mark.parametrize("pc", [0, 1])
def test_donor_results_match_separate_archived_python(tmp_path, unit, phased, dispersion_scope, pc):
    base = bulk_counts().drop(columns="gene")
    counts = pd.concat([base.assign(donor_id="donorA"), base.assign(donor_id="donorB")])
    count_file = tmp_path / "counts.tsv"
    counts.to_csv(count_file, sep="\t", index=False)
    paths = run_per_donor_analysis(
        count_file,
        tmp_path / "results.tsv",
        unit=unit,
        model="single",
        dispersion_scope=dispersion_scope,
        region_col="region" if unit == "feature" else None,
        phased=phased,
        pseudocount=pc,
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t")
    reference_unit = "atac" if unit == "feature" else "snv"
    expected = pd.DataFrame(REFERENCE["bulk"][f"{reference_unit}/{phased}/{pc}"])
    assert set(actual["donor_id"]) == {"donorA", "donorB"}
    for _, donor_results in actual.groupby("donor_id"):
        donor_results = donor_results.copy()
        if unit == "feature":
            donor_results = donor_results.rename(columns={"feature_id": "region"})
        else:
            donor_results["region"] = (
                donor_results["chrom"] + "_" + donor_results["pos"].astype(str)
            )
        comparison = compare_frames(expected, donor_results)
        assert comparison["status"] == "PASS", comparison
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["analysis"]["pvalue_method"] == "chi_square_lrt"


@pytest.mark.rust
def test_donor_default_pseudocount_matches_python(tmp_path):
    counts = bulk_counts().drop(columns="gene").assign(donor_id="donorA")
    count_file = tmp_path / "counts.tsv"
    counts.to_csv(count_file, sep="\t", index=False)
    paths = run_per_donor_analysis(
        count_file,
        tmp_path / "default.tsv",
        unit="feature",
        region_col="region",
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t").rename(columns={"feature_id": "region"})
    comparison = compare_frames(pd.DataFrame(REFERENCE["bulk"]["atac/False/1"]), actual)
    assert comparison["status"] == "PASS", comparison
