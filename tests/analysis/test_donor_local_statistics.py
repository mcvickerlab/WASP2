from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize_scalar
from scipy.special import expit, logsumexp
from scipy.stats import betabinom, chi2, false_discovery_control

from analysis import run_analysis_per_donor as donor_analysis

RHO_EPSILON = 1e-10


def _statistical_counts() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    counts = {
        "donor_a": [(18, 2), (22, 8), (32, 8), (14, 6), (24, 26), (11, 9)],
        "donor_b": [(16, 4), (18, 12), (29, 11), (8, 12), (27, 23), (9, 11)],
    }
    for donor_id, donor_counts in counts.items():
        for index, (ref_count, alt_count) in enumerate(donor_counts, start=1):
            rows.append(
                {
                    "donor_id": donor_id,
                    "chrom": "chr1",
                    "pos": index * 100,
                    "ref": "A",
                    "alt": "G",
                    "GT": "0|1" if index % 2 else "1|0",
                    "region": f"feature_{(index + 1) // 2}",
                    "ref_count": ref_count,
                    "alt_count": alt_count,
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture
def rust_backend(monkeypatch: pytest.MonkeyPatch):
    wasp2_rust = pytest.importorskip("wasp2_rust")
    required = ["fit_imbalance_dispersion", "analyze_imbalance_run", "analyze_imbalance"]
    if not all(hasattr(wasp2_rust, name) for name in required):
        pytest.skip("worktree Rust extension needs rebuilding")
    monkeypatch.setattr(
        donor_analysis,
        "rust_fit_imbalance_dispersion",
        wasp2_rust.fit_imbalance_dispersion,
    )
    monkeypatch.setattr(
        donor_analysis,
        "rust_analyze_imbalance_run",
        wasp2_rust.analyze_imbalance_run,
    )
    return wasp2_rust


def _rho_at_depth(dispersion: pd.Series, total: int) -> float:
    if dispersion["model"] == "single":
        rho = float(dispersion["rho"])
    else:
        center = float(dispersion["linear_depth_center"])
        scale = float(dispersion["linear_depth_scale"])
        assert np.isfinite(center)
        assert np.isfinite(scale) and scale > 0
        standardized_depth = (total - center) / scale
        linear = float(dispersion["linear_d1"]) + standardized_depth * float(
            dispersion["linear_d2"]
        )
        rho = float(expit(linear))
    return float(np.clip(rho, RHO_EPSILON, 1.0 - RHO_EPSILON))


def _logpmf(probability: float, rho: float, ref_count: int, total: int) -> float:
    if probability <= 0.0 or probability >= 1.0:
        return -np.inf
    alpha = probability * (1.0 - rho) / rho
    beta = (1.0 - probability) * (1.0 - rho) / rho
    return float(betabinom.logpmf(ref_count, total, alpha, beta))


def _exact_two_sided_pvalue(ref_count: int, total: int, rho: float) -> float:
    outcomes = np.arange(total + 1)
    alpha = 0.5 * (1.0 - rho) / rho
    log_probabilities = betabinom.logpmf(outcomes, total, alpha, alpha)
    observed_deviation = abs(2 * ref_count - total)
    tail = np.abs(2 * outcomes - total) >= observed_deviation
    return float(np.exp(logsumexp(log_probabilities[tail]) - logsumexp(log_probabilities)))


def _single_snv_oracle(row: pd.Series, dispersion: pd.Series) -> dict[str, float]:
    ref_count = int(row["ref_count"])
    alt_count = int(row["alt_count"])
    total = ref_count + alt_count
    rho = _rho_at_depth(dispersion, total)
    null_ll = _logpmf(0.5, rho, ref_count, total)
    fit = minimize_scalar(
        lambda probability: -_logpmf(probability, rho, ref_count, total),
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    assert fit.success
    alt_ll = -float(fit.fun)
    lrt = max(0.0, 2.0 * (alt_ll - null_ll))
    return {
        "null_ll": null_ll,
        "alt_ll": alt_ll,
        "mu": float(fit.x),
        "lrt": lrt,
        "pval": _exact_two_sided_pvalue(ref_count, total, rho),
    }


def _unphased_feature_oracle(frame: pd.DataFrame, dispersion: pd.Series) -> dict[str, float]:
    observations = []
    for row in frame.sort_values(["chrom", "pos", "ref_count", "alt_count"]).itertuples():
        total = int(row.ref_count + row.alt_count)
        observations.append((int(row.ref_count), total, _rho_at_depth(dispersion, total)))
    null_ll = sum(_logpmf(0.5, rho, ref_count, total) for ref_count, total, rho in observations)

    def negative_log_likelihood(probability: float) -> float:
        log_likelihood = 0.0
        for ref_count, total, rho in observations:
            orientations = [
                _logpmf(probability, rho, ref_count, total),
                _logpmf(1.0 - probability, rho, ref_count, total),
            ]
            log_likelihood += float(logsumexp(orientations) - np.log(2.0))
        return -log_likelihood

    fit = minimize_scalar(
        negative_log_likelihood,
        bounds=(0.5, 1.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    assert fit.success
    alt_ll = -float(fit.fun)
    lrt = max(0.0, 2.0 * (alt_ll - null_ll))
    return {
        "null_ll": null_ll,
        "alt_ll": alt_ll,
        "mu": float(fit.x),
        "lrt": lrt,
        "pval": float(chi2.sf(lrt, 1)),
    }


def _oracle_results(
    counts: pd.DataFrame,
    dispersion: pd.DataFrame,
    *,
    unit: str,
) -> pd.DataFrame:
    rows = []
    for donor_id, donor_counts in counts.groupby("donor_id", sort=True):
        donor_dispersion = dispersion.loc[dispersion["donor_id"] == donor_id].iloc[0]
        if unit == "snv":
            for row in donor_counts.sort_values(["chrom", "pos"]).to_dict("records"):
                series = pd.Series(row)
                rows.append(
                    {
                        "donor_id": donor_id,
                        "identity": (f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"),
                        **_single_snv_oracle(series, donor_dispersion),
                    }
                )
        else:
            for feature_id, feature in donor_counts.groupby("region", sort=True):
                rows.append(
                    {
                        "donor_id": donor_id,
                        "identity": feature_id,
                        **_unphased_feature_oracle(feature, donor_dispersion),
                    }
                )
    oracle = pd.DataFrame(rows)
    oracle["fdr_pval"] = oracle.groupby("donor_id", sort=False)["pval"].transform(
        lambda values: false_discovery_control(values.to_numpy(), method="bh")
    )
    return oracle.sort_values(["donor_id", "identity"]).reset_index(drop=True)


@pytest.mark.rust
@pytest.mark.parametrize("unit,region_col", [("snv", None), ("feature", "region")])
@pytest.mark.parametrize("dispersion_scope", ["global", "per-donor"])
@pytest.mark.parametrize("model", ["single", "linear"])
def test_model_scope_unit_matrix_matches_independent_scipy_oracle(
    tmp_path: Path,
    rust_backend: Any,
    unit: str,
    region_col: str | None,
    dispersion_scope: str,
    model: str,
) -> None:
    del rust_backend
    counts = _statistical_counts()
    count_path = tmp_path / "counts.tsv"
    counts.to_csv(count_path, sep="\t", index=False)
    paths = donor_analysis.run_per_donor_analysis(
        count_path,
        tmp_path / f"{unit}-{model}-{dispersion_scope}.tsv",
        unit=unit,
        model=model,
        dispersion_scope=dispersion_scope,
        region_col=region_col,
        phased=False,
        min_count=0,
        pseudocount=0,
        min_donor_observations=4,
    )

    observed = pd.read_csv(paths["results"], sep="\t")
    identity_column = "snv_id" if unit == "snv" else "feature_id"
    observed = observed.rename(columns={identity_column: "identity"}).sort_values(
        ["donor_id", "identity"]
    )
    observed = observed.reset_index(drop=True)
    dispersion = pd.read_csv(paths["dispersion"], sep="\t")
    oracle = _oracle_results(counts, dispersion, unit=unit)

    assert observed[["donor_id", "identity"]].equals(oracle[["donor_id", "identity"]])
    for column in ["null_ll", "alt_ll", "mu", "lrt", "pval", "fdr_pval"]:
        np.testing.assert_allclose(
            observed[column],
            oracle[column],
            rtol=3e-5,
            atol=3e-5,
            err_msg=f"{model}/{dispersion_scope}/{unit}: {column}",
        )
    assert observed["ref_count"].sum() == counts["ref_count"].sum()
    assert observed["alt_count"].sum() == counts["alt_count"].sum()
    if dispersion_scope == "global":
        assert dispersion["fit_id"].nunique() == 1
        parameters = (
            ["rho"]
            if model == "single"
            else [
                "linear_d1",
                "linear_d2",
                "linear_depth_center",
                "linear_depth_scale",
            ]
        )
        assert all(dispersion[column].nunique() == 1 for column in parameters)
    else:
        assert dispersion["fit_id"].nunique() == 2


def _run_legacy(
    rust_backend: Any,
    path: Path,
    *,
    method: str,
    phased: bool,
    region_col: str | None,
) -> pd.DataFrame:
    results = rust_backend.analyze_imbalance(
        str(path),
        min_count=0,
        pseudocount=0,
        method=method,
        phased=phased,
        region_col=region_col,
    )
    return pd.DataFrame(results).sort_values("region").reset_index(drop=True)


@pytest.mark.rust
@pytest.mark.parametrize("method", ["single", "linear"])
def test_unphased_feature_results_are_row_permutation_invariant(
    tmp_path: Path,
    rust_backend: Any,
    method: str,
) -> None:
    counts = (
        _statistical_counts()
        .loc[lambda frame: frame["donor_id"] == "donor_a"]
        .drop(columns="donor_id")
    )
    forward = tmp_path / "forward.tsv"
    shuffled = tmp_path / "shuffled.tsv"
    counts.to_csv(forward, sep="\t", index=False)
    counts.sample(frac=1.0, random_state=1729).to_csv(shuffled, sep="\t", index=False)

    first = _run_legacy(rust_backend, forward, method=method, phased=False, region_col="region")
    second = _run_legacy(rust_backend, shuffled, method=method, phased=False, region_col="region")
    assert first[["region", "ref_count", "alt_count", "N", "snp_count"]].equals(
        second[["region", "ref_count", "alt_count", "N", "snp_count"]]
    )
    for column in ["null_ll", "alt_ll", "mu", "lrt", "pval", "fdr_pval"]:
        np.testing.assert_allclose(first[column], second[column], rtol=0, atol=1e-12)


@pytest.mark.rust
@pytest.mark.parametrize("method", ["single", "linear"])
def test_allele_swap_preserves_tests_and_mirrors_single_snp_mu(
    tmp_path: Path,
    rust_backend: Any,
    method: str,
) -> None:
    counts = (
        _statistical_counts()
        .loc[lambda frame: frame["donor_id"] == "donor_a"]
        .drop(columns=["donor_id", "region"])
    )
    original_path = tmp_path / "original.tsv"
    swapped_path = tmp_path / "swapped.tsv"
    counts.to_csv(original_path, sep="\t", index=False)
    swapped = counts.copy()
    swapped[["ref", "alt"]] = swapped[["alt", "ref"]].to_numpy()
    swapped[["ref_count", "alt_count"]] = swapped[["alt_count", "ref_count"]].to_numpy()
    swapped["GT"] = swapped["GT"].map({"0|1": "1|0", "1|0": "0|1"})
    swapped.to_csv(swapped_path, sep="\t", index=False)

    original = _run_legacy(
        rust_backend, original_path, method=method, phased=False, region_col=None
    )
    swapped_results = _run_legacy(
        rust_backend, swapped_path, method=method, phased=False, region_col=None
    )
    for column in ["null_ll", "alt_ll", "lrt", "pval", "fdr_pval"]:
        np.testing.assert_allclose(original[column], swapped_results[column], rtol=0, atol=2e-10)
    np.testing.assert_allclose(
        original["mu"] + swapped_results["mu"], np.ones(len(original)), rtol=0, atol=2e-5
    )


@pytest.mark.rust
@pytest.mark.parametrize("method", ["single", "linear"])
def test_single_snp_results_are_phase_invariant(
    tmp_path: Path,
    rust_backend: Any,
    method: str,
) -> None:
    counts = (
        _statistical_counts()
        .loc[lambda frame: frame["donor_id"] == "donor_a"]
        .drop(columns=["donor_id", "region"])
    )
    path = tmp_path / "counts.tsv"
    counts.to_csv(path, sep="\t", index=False)
    unphased = _run_legacy(rust_backend, path, method=method, phased=False, region_col=None)
    phased = _run_legacy(rust_backend, path, method=method, phased=True, region_col=None)

    assert unphased[["region", "ref_count", "alt_count", "N", "snp_count"]].equals(
        phased[["region", "ref_count", "alt_count", "N", "snp_count"]]
    )
    for column in ["null_ll", "alt_ll", "mu", "lrt", "pval", "fdr_pval"]:
        np.testing.assert_allclose(unphased[column], phased[column], rtol=0, atol=1e-12)


@pytest.mark.rust
@pytest.mark.parametrize("method", ["per_donor", "per-donor", "single-global", "mystery"])
def test_rust_rejects_unknown_and_retired_pooled_methods(
    tmp_path: Path,
    rust_backend: Any,
    method: str,
) -> None:
    path = tmp_path / "counts.tsv"
    _statistical_counts().drop(columns="donor_id").to_csv(path, sep="\t", index=False)
    with pytest.raises(RuntimeError, match="removed|Unknown method"):
        rust_backend.analyze_imbalance(str(path), method=method)


@pytest.mark.rust
def test_rust_has_no_cohort_shared_entry_point(rust_backend: Any) -> None:
    assert not hasattr(rust_backend, "analyze_cohort_snvs")
