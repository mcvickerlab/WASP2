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


@pytest.mark.rust
@pytest.mark.parametrize("dispersion_scope", ["global", "per-donor"])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
def test_overlapping_donor_features_match_python(tmp_path, dispersion_scope, phased, pc):
    """Preserve Python observation weights with heterogeneous donors and overlaps."""
    import contextlib
    import io

    import numpy as np
    from scipy.optimize import OptimizeResult, minimize_scalar
    from scipy.stats import betabinom
    from scripts.compare_archive_analysis import load_archive

    first = bulk_counts().drop(columns="gene")
    second = first.copy()
    second["ref_count"] = [11, 12, 21, 20, 17, 22, 8, 11, 19, 21, 9, 12]
    second["alt_count"] = [9, 8, 19, 20, 23, 18, 12, 9, 21, 19, 11, 8]
    # An identical observation in both donors must survive the global fit.
    second.loc[2, ["ref_count", "alt_count"]] = first.loc[2, ["ref_count", "alt_count"]]
    donors = {}
    for donor, base in [("donorA", first), ("donorB", second)]:
        extra = base.iloc[:2].copy()
        extra["region"] = "overlap"
        # Reverse order to exercise the original first-SNP anchor as well.
        donors[donor] = pd.concat([base, extra], ignore_index=True).iloc[::-1].copy()
    retained = pd.concat([base.assign(donor_id=donor) for donor, base in donors.items()])
    counts = pd.concat([retained, retained.iloc[:1]], ignore_index=True)
    count_file = tmp_path / "counts.tsv"
    counts.to_csv(count_file, sep="\t", index=False)
    paths = run_per_donor_analysis(
        count_file,
        tmp_path / "results.tsv",
        unit="feature",
        region_col="region",
        dispersion_scope=dispersion_scope,
        phased=phased,
        pseudocount=pc,
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t").rename(columns={"feature_id": "region"})
    fits = pd.read_csv(paths["dispersion"], sep="\t").set_index("donor_id")
    qc = pd.read_csv(paths["qc"], sep="\t").set_index("donor_id")
    archive = load_archive()["as_analysis"]
    original_optimizer = archive.minimize_scalar
    global_rho = None
    if dispersion_scope == "global":
        ref = retained["ref_count"].to_numpy() + pc
        depth = retained["ref_count"].to_numpy() + retained["alt_count"].to_numpy() + 2 * pc

        def balanced_null(rho):
            shape = 0.5 * (1 - rho) / rho
            return -np.sum(betabinom.logpmf(ref, depth, shape, shape))

        global_rho = minimize_scalar(balanced_null, method="bounded", bounds=(0, 1)).x

        def shared_null_optimizer(fun, *args, **kwargs):
            # Only replace the balanced-null dispersion fit. Effect optimizers
            # still execute the archived Python functions separately per donor.
            if len(kwargs.get("args", ())) == 2:
                return OptimizeResult(x=global_rho)
            return original_optimizer(fun, *args, **kwargs)

        archive.minimize_scalar = shared_null_optimizer
    for donor, base in donors.items():
        with contextlib.redirect_stdout(io.StringIO()):
            expected = archive.get_imbalance(
                base.copy(),
                region_col="region",
                phased=phased,
                pseudocount=pc,
            )
        comparison = compare_frames(expected, actual.loc[actual["donor_id"] == donor])
        assert comparison["status"] == "PASS", comparison
        assert fits.loc[donor, "fit_observations"] == (28 if dispersion_scope == "global" else 14)
        assert fits.loc[donor, "donor_eligible_snv_observations"] == 12
        assert qc.loc[donor, "eligible_snv_observations"] == 12
        assert qc.loc[donor, "eligible_test_rows"] == 14
        if global_rho is not None:
            np.testing.assert_allclose(fits.loc[donor, "rho"], global_rho, rtol=5e-5, atol=5e-5)


@pytest.mark.rust
@pytest.mark.parametrize("scope", ["global", "per-donor"])
@pytest.mark.parametrize("pc", [0, 1])
@pytest.mark.parametrize(
    "case,reason",
    [
        ("missing", "missing_gt"),
        ("same_forward", "identical_gt"),
        ("same_reverse", "identical_gt"),
        ("alleles", "no_numeric_phased_gt"),
        ("unphased", "no_numeric_phased_gt"),
        ("filtered", "identical_gt"),
        ("informative", None),
    ],
)
def test_donor_phasing_matches_archive_fallback(tmp_path, scope, pc, case, reason):
    import contextlib
    import io

    from scripts.compare_archive_analysis import load_archive

    base = bulk_counts().drop(columns="gene")
    if case == "missing":
        base = base.drop(columns="GT")
    elif case in {"same_forward", "same_reverse"}:
        base["GT"] = "0|1" if case == "same_forward" else "1|0"
    elif case == "alleles":
        base["GT"] = base["GT"].replace({"0|1": "A|G", "1|0": "G|A"})
    elif case == "unphased":
        base["GT"] = base["GT"].replace({"0|1": "0/1", "1|0": "1/0"})
    elif case == "filtered":
        base.loc[base["GT"] == "1|0", ["ref_count", "alt_count"]] = 0
    file = tmp_path / "counts.tsv"
    base.assign(donor_id="donorA").to_csv(file, sep="\t", index=False)
    with contextlib.redirect_stdout(io.StringIO()):
        expected = load_archive()["as_analysis"].get_imbalance(
            base.copy(),
            region_col="region",
            phased=True,
            pseudocount=pc,
        )
    paths = run_per_donor_analysis(
        file,
        tmp_path / "results.tsv",
        unit="feature",
        region_col="region",
        phased=True,
        pseudocount=pc,
        dispersion_scope=scope,
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t").rename(columns={"feature_id": "region"})
    comparison = compare_frames(expected, actual)
    assert comparison["status"] == "PASS", comparison
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["analysis"]["phased"] is True
    assert provenance["donors"]["phasing"] == [
        {
            "donor_id": "donorA",
            "requested_phased": True,
            "effective_phased": reason is None,
            "fallback_reason": reason,
        }
    ]


@pytest.mark.rust
def test_donor_phasing_decisions_are_independent(tmp_path):
    base = bulk_counts().drop(columns="gene")
    first = base.assign(donor_id="informative")
    second = base.assign(donor_id="fallback", GT="0|1")
    file = tmp_path / "counts.tsv"
    pd.concat([first, second]).to_csv(file, sep="\t", index=False)
    paths = run_per_donor_analysis(
        file,
        tmp_path / "results.tsv",
        unit="feature",
        region_col="region",
        phased=True,
        dispersion_scope="global",
        min_donor_observations=4,
    )
    provenance = json.loads(paths["provenance"].read_text())
    decisions = {row["donor_id"]: row for row in provenance["donors"]["phasing"]}
    assert decisions["informative"]["effective_phased"] is True
    assert decisions["fallback"]["effective_phased"] is False
    assert decisions["fallback"]["fallback_reason"] == "identical_gt"


@pytest.mark.rust
@pytest.mark.parametrize("unit,phased", [("snv", False), ("feature", False), ("feature", True)])
@pytest.mark.parametrize("pc", [0, 1])
@pytest.mark.parametrize("scope", ["global", "per-donor"])
def test_distinct_repeated_snp_observations_match_archive(tmp_path, unit, phased, pc, scope):
    import contextlib
    import io

    from scripts.compare_archive_analysis import load_archive

    base = bulk_counts().drop(columns="gene")
    changed = base.iloc[:1].copy()
    changed["ref_count"] += 3
    opposite = base.iloc[:1].copy()
    opposite["GT"] = "1|0"
    parts = [base, changed, opposite, changed.copy()]
    if unit == "feature":
        overlapping = changed.copy()
        overlapping["region"] = "overlap"
        parts.append(overlapping)
    base = pd.concat(parts, ignore_index=True)
    if unit == "snv":
        base = base.drop(columns="region")
    counts = pd.concat([base.assign(donor_id=donor) for donor in ["A", "B"]])
    file = tmp_path / "counts.tsv"
    counts.to_csv(file, sep="\t", index=False)
    with contextlib.redirect_stdout(io.StringIO()):
        expected = load_archive()["as_analysis"].get_imbalance(
            base.copy(),
            region_col="region" if unit == "feature" else None,
            phased=phased,
            pseudocount=pc,
        )
    if unit == "snv":
        expected = expected.rename(columns={"variant": "region"})
    paths = run_per_donor_analysis(
        file,
        tmp_path / "results.tsv",
        unit=unit,
        region_col="region" if unit == "feature" else None,
        phased=phased,
        pseudocount=pc,
        dispersion_scope=scope,
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t")
    for _, observed in actual.groupby("donor_id"):
        observed = observed.copy()
        if unit == "feature":
            observed = observed.rename(columns={"feature_id": "region"})
        else:
            observed["region"] = observed["chrom"] + "_" + observed["pos"].astype(str)
        comparison = compare_frames(expected, observed)
        assert comparison["status"] == "PASS", comparison


@pytest.mark.rust
def test_fallback_duplicate_gt_is_ignored_after_phase_decision(tmp_path):
    import contextlib
    import io

    from scripts.compare_archive_analysis import load_archive

    base = bulk_counts().drop(columns="gene").assign(GT="A|G")
    duplicate = base.iloc[:1].copy()
    duplicate["GT"] = "G|A"
    base = pd.concat([base, duplicate], ignore_index=True)
    file = tmp_path / "counts.tsv"
    base.assign(donor_id="A").to_csv(file, sep="\t", index=False)
    with contextlib.redirect_stdout(io.StringIO()):
        expected = load_archive()["as_analysis"].get_imbalance(
            base.copy(),
            region_col="region",
            phased=True,
        )
    paths = run_per_donor_analysis(
        file,
        tmp_path / "results.tsv",
        unit="feature",
        region_col="region",
        phased=True,
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t").rename(columns={"feature_id": "region"})
    assert compare_frames(expected, actual)["status"] == "PASS"
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["donors"]["phasing"][0]["fallback_reason"] == "no_numeric_phased_gt"
    assert pd.read_csv(paths["dispersion"], sep="\t")["fit_observations"].item() == 12
