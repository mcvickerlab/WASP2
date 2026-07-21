from __future__ import annotations

import hashlib
import importlib.util
import json
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control
from typer.testing import CliRunner

from analysis import __main__ as analysis_cli
from analysis import run_analysis
from analysis import run_analysis_per_donor as donor_analysis


def _count_frame(*, overlap: bool = False, exact_duplicate: bool = False) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    donor_counts = {
        "donor_a": [(18, 2), (10, 10), (14, 6), (7, 13)],
        "donor_b": [(16, 4), (9, 11), (13, 7), (8, 12)],
    }
    for donor_id, counts in donor_counts.items():
        for index, (ref_count, alt_count) in enumerate(counts, start=1):
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
    if overlap:
        for donor_id in donor_counts:
            repeated = next(row.copy() for row in rows if row["donor_id"] == donor_id)
            repeated["region"] = "feature_overlap"
            rows.append(repeated)
    if exact_duplicate:
        rows.append(rows[0].copy())
    return pd.DataFrame(rows)


def _write_counts(
    path: Path,
    *,
    donor_column: str = "donor_id",
    overlap: bool = False,
    exact_duplicate: bool = False,
) -> pd.DataFrame:
    frame = _count_frame(overlap=overlap, exact_duplicate=exact_duplicate)
    if donor_column != "donor_id":
        frame = frame.rename(columns={"donor_id": donor_column})
    frame.to_csv(path, sep="\t", index=False)
    return frame


def _default_pvalues(frame: pd.DataFrame, groups: list[tuple[str, pd.DataFrame]]) -> list[float]:
    del frame
    return [min(0.95, 0.01 * (index + 1)) for index in range(len(groups))]


def _install_backend(
    monkeypatch: pytest.MonkeyPatch,
    *,
    pvalues: Callable[[pd.DataFrame, list[tuple[str, pd.DataFrame]]], list[float]] = (
        _default_pvalues
    ),
    scalar_rho: float = float.fromhex("0x1.23456789abcdep-5"),
    linear_d1: float = float.fromhex("-0x1.3456789abcdep+1"),
    linear_d2: float = float.fromhex("-0x1.23456789abcdep-6"),
) -> dict[str, list[dict[str, Any]]]:
    calls: dict[str, list[dict[str, Any]]] = {"fit": [], "analyze": []}

    def fake_fit(path: str, *, min_count: int, pseudocount: int, method: str):
        frame = pd.read_csv(path, sep="\t")
        donors = tuple(frame["donor_id"].drop_duplicates())
        fit_id = f"{method}:{','.join(donors)}"
        payload = {
            "fit_id": fit_id,
            "method": method,
            "rho": scalar_rho if method == "single" else None,
            "linear_d1": linear_d1 if method == "linear" else None,
            "linear_d2": linear_d2 if method == "linear" else None,
            "n_observations": len(frame),
        }
        calls["fit"].append(
            {
                "frame": frame.copy(),
                "min_count": min_count,
                "pseudocount": pseudocount,
                "method": method,
                "payload": payload,
            }
        )
        return payload

    def fake_analyze(
        path: str,
        *,
        min_count: int,
        pseudocount: int,
        phased: bool,
        region_col: str | None,
        method: str,
        dispersion_fit: dict[str, Any],
    ):
        frame = pd.read_csv(path, sep="\t")
        if region_col is None:
            groups = [
                (f"{row.chrom}_{row.pos}", frame.iloc[[index]])
                for index, row in enumerate(frame.itertuples(index=False))
            ]
        else:
            groups = [
                (str(region), group.copy())
                for region, group in frame.groupby(region_col, sort=False)
            ]
        raw_pvalues = pvalues(frame, groups)
        adjusted = false_discovery_control(raw_pvalues, method="bh")
        results = []
        for (region, group), pval, fdr_pval in zip(groups, raw_pvalues, adjusted):
            ref_count = int(group["ref_count"].sum())
            alt_count = int(group["alt_count"].sum())
            results.append(
                {
                    "region": region,
                    "ref_count": ref_count,
                    "alt_count": alt_count,
                    "N": ref_count + alt_count,
                    "snp_count": len(group),
                    "null_ll": -3.0,
                    "alt_ll": -2.0,
                    "mu": 0.6,
                    "lrt": 2.0,
                    "pval": float(pval),
                    "fdr_pval": float(fdr_pval),
                }
            )
        calls["analyze"].append(
            {
                "frame": frame.copy(),
                "min_count": min_count,
                "pseudocount": pseudocount,
                "phased": phased,
                "region_col": region_col,
                "method": method,
                "dispersion_fit": dispersion_fit,
            }
        )
        return {
            "results": results,
            "n_observations": len(frame),
            "requested_phased": phased,
            "effective_phased": phased,
            "dispersion_fit": dispersion_fit,
        }

    monkeypatch.setattr(donor_analysis, "rust_fit_imbalance_dispersion", fake_fit)
    monkeypatch.setattr(donor_analysis, "rust_analyze_imbalance_run", fake_analyze)
    return calls


@pytest.mark.unit
@pytest.mark.parametrize("unit,region_col", [("snv", None), ("feature", "region")])
@pytest.mark.parametrize("dispersion_scope", ["global", "per-donor"])
@pytest.mark.parametrize("model", ["single", "linear"])
def test_model_scope_unit_matrix(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    unit: str,
    region_col: str | None,
    dispersion_scope: str,
    model: str,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    calls = _install_backend(monkeypatch)

    paths = donor_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "results.tsv",
        unit=unit,
        model=model,
        dispersion_scope=dispersion_scope,
        region_col=region_col,
        phased=False,
        min_count=0,
        pseudocount=0,
        min_donor_observations=1,
    )

    assert len(calls["fit"]) == (1 if dispersion_scope == "global" else 2)
    assert len(calls["analyze"]) == 2
    assert {call["method"] for call in calls["fit"]} == {model}
    assert {call["method"] for call in calls["analyze"]} == {model}
    assert {call["region_col"] for call in calls["analyze"]} == {region_col}
    assert all("donor_id" not in call["frame"] for call in calls["analyze"])
    assert all(call["pseudocount"] == 0 for call in calls["fit"] + calls["analyze"])

    fit_ids = [call["dispersion_fit"]["fit_id"] for call in calls["analyze"]]
    assert (len(set(fit_ids)) == 1) is (dispersion_scope == "global")
    results = pd.read_csv(paths["results"], sep="\t")
    identity = "snv_id" if unit == "snv" else "feature_id"
    assert results.groupby("donor_id")[identity].count().to_dict() == {
        "donor_a": 4 if unit == "snv" else 2,
        "donor_b": 4 if unit == "snv" else 2,
    }
    dispersion = pd.read_csv(paths["dispersion"], sep="\t")
    assert dispersion["model"].eq(model).all()
    assert dispersion["dispersion_scope"].eq(dispersion_scope).all()


@pytest.mark.unit
@pytest.mark.parametrize("model", ["single", "linear"])
def test_global_fit_payload_is_reused_exactly(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    model: str,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    calls = _install_backend(monkeypatch)

    donor_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "results.tsv",
        unit="snv",
        model=model,
        dispersion_scope="global",
        min_count=0,
        pseudocount=0,
        min_donor_observations=1,
    )

    fitted = calls["fit"][0]["payload"]
    reused = [call["dispersion_fit"] for call in calls["analyze"]]
    assert reused[0] is reused[1]
    assert reused[0] == fitted
    assert reused[0]["fit_id"] == fitted["fit_id"]
    parameter_names = ["rho"] if model == "single" else ["linear_d1", "linear_d2"]
    for parameter in parameter_names:
        assert float(reused[0][parameter]).hex() == float(fitted[parameter]).hex()
        assert float(reused[1][parameter]).hex() == float(fitted[parameter]).hex()


@pytest.mark.unit
def test_bh_correction_and_likelihood_calls_are_donor_local(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    expected = {
        "donor_a": np.array([0.001, 0.04, 0.2, 0.8]),
        "donor_b": np.array([0.03, 0.031, 0.032, 0.9]),
    }

    def donor_pvalues(frame: pd.DataFrame, groups: list[tuple[str, pd.DataFrame]]) -> list[float]:
        donor_id = "donor_a" if int(frame.iloc[0]["ref_count"]) == 18 else "donor_b"
        assert len(groups) == len(expected[donor_id])
        return expected[donor_id].tolist()

    calls = _install_backend(monkeypatch, pvalues=donor_pvalues)
    paths = donor_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "results.tsv",
        unit="snv",
        model="single",
        dispersion_scope="per-donor",
        min_count=0,
        pseudocount=0,
        min_donor_observations=1,
    )

    assert len(calls["analyze"]) == 2
    assert all(len(call["frame"]) == 4 for call in calls["analyze"])
    results = pd.read_csv(paths["results"], sep="\t")
    observed_qvalues = []
    for donor_id, raw_pvalues in expected.items():
        donor_results = results.loc[results["donor_id"] == donor_id].sort_values("snv_id")
        np.testing.assert_allclose(donor_results["pval"], raw_pvalues, rtol=0, atol=0)
        expected_qvalues = false_discovery_control(raw_pvalues, method="bh")
        np.testing.assert_allclose(donor_results["fdr_pval"], expected_qvalues, rtol=0, atol=1e-15)
        observed_qvalues.extend(donor_results["fdr_pval"].tolist())
    pooled_qvalues = false_discovery_control(np.concatenate(list(expected.values())), method="bh")
    assert not np.allclose(observed_qvalues, pooled_qvalues)


@pytest.mark.unit
@pytest.mark.parametrize("donor_column", ["donor_id", "sample"])
def test_canonical_and_legacy_donor_columns_normalize_to_donor_id(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    donor_column: str,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts, donor_column=donor_column)
    _install_backend(monkeypatch)

    paths = donor_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "results.tsv",
        unit="snv",
        min_count=0,
        pseudocount=0,
        min_donor_observations=1,
    )

    results = pd.read_csv(paths["results"], sep="\t")
    assert "donor_id" in results
    assert "sample" not in results
    assert sorted(results["donor_id"].unique()) == ["donor_a", "donor_b"]
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["inputs"]["qc"]["legacy_sample_normalized"] is (donor_column == "sample")


@pytest.mark.unit
def test_both_donor_columns_fail_closed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    counts = tmp_path / "counts.tsv"
    frame = _write_counts(counts)
    frame["sample"] = frame["donor_id"]
    frame.to_csv(counts, sep="\t", index=False)
    _install_backend(monkeypatch)

    with pytest.raises(ValueError, match="both donor_id and legacy sample"):
        donor_analysis.run_per_donor_analysis(
            counts,
            tmp_path / "results.tsv",
            unit="snv",
            min_donor_observations=1,
        )


@pytest.mark.unit
def test_overlapping_feature_membership_is_preserved_but_fit_once_per_snv(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts, overlap=True, exact_duplicate=True)
    calls = _install_backend(monkeypatch)

    paths = donor_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "results.tsv",
        unit="feature",
        model="single",
        dispersion_scope="global",
        region_col="region",
        min_count=0,
        pseudocount=0,
        min_donor_observations=1,
    )

    assert len(calls["fit"]) == 1
    assert len(calls["fit"][0]["frame"]) == 8
    assert [len(call["frame"]) for call in calls["analyze"]] == [5, 5]
    results = pd.read_csv(paths["results"], sep="\t")
    assert results.groupby("donor_id")["feature_id"].apply(set).to_dict() == {
        "donor_a": {"feature_1", "feature_2", "feature_overlap"},
        "donor_b": {"feature_1", "feature_2", "feature_overlap"},
    }
    provenance = json.loads(paths["provenance"].read_text())
    input_qc = provenance["inputs"]["qc"]
    assert input_qc["raw_rows"] == 11
    assert input_qc["validated_rows"] == 10
    assert input_qc["exact_duplicate_rows_removed"] == 1
    assert input_qc["unique_donor_snv_rows"] == 8


@pytest.mark.unit
@pytest.mark.parametrize("conflict_column,replacement", [("alt", "T"), ("ref_count", 17)])
def test_conflicting_repeated_feature_rows_fail_closed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    conflict_column: str,
    replacement: Any,
) -> None:
    counts = tmp_path / "counts.tsv"
    frame = _count_frame(overlap=True)
    overlap_index = frame.index[
        (frame["donor_id"] == "donor_a") & (frame["region"] == "feature_overlap")
    ][0]
    frame.loc[overlap_index, conflict_column] = replacement
    frame.to_csv(counts, sep="\t", index=False)
    calls = _install_backend(monkeypatch)

    with pytest.raises(ValueError, match="Conflicting"):
        donor_analysis.run_per_donor_analysis(
            counts,
            tmp_path / "results.tsv",
            unit="feature",
            region_col="region",
            min_count=0,
            min_donor_observations=1,
        )
    assert calls == {"fit": [], "analyze": []}


@pytest.mark.unit
def test_phased_snv_rejected_before_input_io(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _install_backend(monkeypatch)
    with pytest.raises(ValueError, match="SNV analysis is unphased"):
        donor_analysis.run_per_donor_analysis(
            "does-not-exist.tsv",
            "results.tsv",
            unit="snv",
            phased=True,
        )
    assert calls == {"fit": [], "analyze": []}


@pytest.mark.unit
@pytest.mark.parametrize("bad_model", ["per_donor", "per-donor", "single-global", "mystery"])
def test_unknown_and_legacy_pooled_models_fail_before_input_io(
    monkeypatch: pytest.MonkeyPatch,
    bad_model: str,
) -> None:
    calls = _install_backend(monkeypatch)
    with pytest.raises(ValueError, match="single.*linear"):
        donor_analysis.run_per_donor_analysis(
            "does-not-exist.tsv",
            "results.tsv",
            unit="snv",
            model=bad_model,
        )
    assert calls == {"fit": [], "analyze": []}


@pytest.mark.unit
def test_per_donor_route_defaults_to_raw_counts_without_pseudocount(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, Any] = {}

    def fake_run(count_file: str | Path, out_file: str | Path, **kwargs: Any):
        observed.update({"count_file": count_file, "out_file": out_file, **kwargs})
        return {"results": Path(out_file)}

    monkeypatch.setattr(run_analysis, "run_per_donor_analysis", fake_run)
    run_analysis.run_ai_analysis(
        "counts.tsv",
        out_file="results.tsv",
        scope="per-donor",
        unit="feature",
        region_col="region",
        model="linear",
        dispersion_scope="global",
    )

    assert observed["unit"] == "feature"
    assert observed["model"] == "linear"
    assert observed["dispersion_scope"] == "global"
    assert observed["pseudocount"] == 0


@pytest.mark.unit
@pytest.mark.parametrize("unit,region_col", [("snv", None), ("feature", "region")])
@pytest.mark.parametrize("dispersion_scope", ["global", "per-donor"])
@pytest.mark.parametrize("model", ["single", "linear"])
def test_run_ai_analysis_forwards_model_scope_unit_matrix(
    monkeypatch: pytest.MonkeyPatch,
    unit: str,
    region_col: str | None,
    dispersion_scope: str,
    model: str,
) -> None:
    observed: dict[str, Any] = {}

    def fake_run(count_file: str | Path, out_file: str | Path, **kwargs: Any):
        observed.update({"count_file": count_file, "out_file": out_file, **kwargs})
        return {"results": Path(out_file)}

    monkeypatch.setattr(run_analysis, "run_per_donor_analysis", fake_run)
    run_analysis.run_ai_analysis(
        "counts.tsv",
        out_file="results.tsv",
        scope="per-donor",
        unit=unit,
        region_col=region_col,
        model=model,
        dispersion_scope=dispersion_scope,
    )

    assert observed["unit"] == unit
    assert observed["region_col"] == region_col
    assert observed["model"] == model
    assert observed["dispersion_scope"] == dispersion_scope
    assert observed["pseudocount"] == 0


@pytest.mark.unit
def test_console_interface_forwards_donor_local_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, Any] = {}

    def fake_run(**kwargs: Any) -> None:
        observed.update(kwargs)

    monkeypatch.setattr(analysis_cli, "run_ai_analysis", fake_run)
    manifest_digest = "a" * 64
    result = CliRunner().invoke(
        analysis_cli.app,
        [
            "find-imbalance",
            "count-bundle",
            "--scope",
            "per-donor",
            "--unit",
            "feature",
            "--model",
            "linear",
            "--dispersion-scope",
            "global",
            "--region-col",
            "region",
            "--pseudocount",
            "0",
            "--min-donor-observations",
            "12",
            "--expected-manifest-sha256",
            manifest_digest,
            "-o",
            "results.tsv",
        ],
    )

    assert result.exit_code == 0, result.output
    assert observed["count_file"] == "count-bundle"
    assert observed["out_file"] == "results.tsv"
    assert observed["scope"] == "per-donor"
    assert observed["unit"] == "feature"
    assert observed["model"] == "linear"
    assert observed["dispersion_scope"] == "global"
    assert observed["pseudocount"] == 0
    assert observed["min_donor_observations"] == 12
    assert observed["expected_manifest_sha256"] == manifest_digest


@pytest.mark.unit
@pytest.mark.parametrize("bad_model", ["per_donor", "per-donor", "single-global", "mystery"])
def test_legacy_global_route_rejects_unknown_models(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    bad_model: str,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    monkeypatch.setattr(run_analysis, "rust_analyze_imbalance", lambda *args, **kwargs: [])

    with pytest.raises(ValueError, match="single.*linear"):
        run_analysis.run_ai_analysis(
            counts,
            out_file=tmp_path / "results.tsv",
            model=bad_model,
            per_variant=True,
        )


@pytest.mark.unit
def test_retired_cohort_shared_surface_is_absent() -> None:
    assert importlib.util.find_spec("analysis.run_analysis_cohort_shared") is None
    assert not hasattr(run_analysis, "run_cohort_shared_snv_analysis")
    with pytest.raises(ValueError, match="scope must be 'per-donor'"):
        run_analysis.run_ai_analysis(
            "does-not-exist.tsv",
            scope="cohort-shared",
            unit="snv",
        )


@pytest.mark.unit
@pytest.mark.parametrize("unit,region_col", [("snv", None), ("feature", "region")])
def test_sidecar_and_provenance_schema(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    unit: str,
    region_col: str | None,
) -> None:
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    _install_backend(monkeypatch)
    paths = donor_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "results.tsv",
        unit=unit,
        model="linear",
        dispersion_scope="global",
        region_col=region_col,
        min_count=0,
        pseudocount=0,
        min_donor_observations=1,
    )

    assert set(paths) == {"results", "dispersion", "qc", "provenance"}
    assert all(path.is_file() for path in paths.values())
    result_prefix = (
        ["donor_id", "snv_id", "chrom", "pos", "ref", "alt"]
        if unit == "snv"
        else ["donor_id", "feature_id"]
    )
    assert pd.read_csv(paths["results"], sep="\t").columns.tolist() == [
        *result_prefix,
        "ref_count",
        "alt_count",
        "N",
        "snp_count",
        "null_ll",
        "alt_ll",
        "mu",
        "lrt",
        "pval",
        "fdr_pval",
        "significant_q05",
        "significant_q10",
    ]
    assert pd.read_csv(paths["dispersion"], sep="\t").columns.tolist() == [
        "donor_id",
        "dispersion_scope",
        "fit_id",
        "model",
        "rho",
        "linear_d1",
        "linear_d2",
        "fit_observations",
        "donor_eligible_snv_observations",
        "result_rows",
    ]
    assert pd.read_csv(paths["qc"], sep="\t").columns.tolist() == [
        "donor_id",
        "input_rows",
        "unique_snv_rows",
        "eligible_snv_observations",
        "eligible_test_rows",
        "analyzed",
        "status",
    ]

    provenance = json.loads(paths["provenance"].read_text())
    assert set(provenance) == {
        "schema_version",
        "created_at_utc",
        "software",
        "analysis",
        "parameters",
        "inputs",
        "donors",
        "outputs",
    }
    assert provenance["schema_version"] == donor_analysis.SCHEMA_VERSION
    datetime.fromisoformat(provenance["created_at_utc"])
    assert provenance["analysis"] == {
        "scope": "per-donor",
        "unit": unit,
        "model": "linear",
        "dispersion_scope": "global",
        "region_col": region_col,
        "phased": False,
        "multiple_testing": "Benjamini-Hochberg independently within each donor",
    }
    assert provenance["parameters"] == {
        "min_count": 0,
        "pseudocount": 0,
        "min_donor_observations": 1,
    }
    assert provenance["donors"] == {"total": 2, "included": 2, "excluded": 0}
    assert set(provenance["outputs"]) == {"results", "dispersion", "qc"}
    for name, record in provenance["outputs"].items():
        assert set(record) == {"path", "sha256", "size_bytes", "mtime_utc"}
        assert record["path"] == paths[name].name
        assert record["sha256"] == hashlib.sha256(paths[name].read_bytes()).hexdigest()
        assert record["size_bytes"] == paths[name].stat().st_size
        datetime.fromisoformat(record["mtime_utc"])
