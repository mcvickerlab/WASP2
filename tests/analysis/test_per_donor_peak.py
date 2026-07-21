import json
from pathlib import Path

import pandas as pd
import pytest

from analysis import run_analysis_per_donor as peak_analysis


def _write_peak_counts(path: Path, *, bad_gt: bool = False) -> None:
    rows = []
    for donor in ["donor_a", "donor_b", "donor_low"]:
        count = 50 if donor != "donor_low" else 2
        for index in range(count):
            ref, alt = ("A", "G")
            gt = f"{ref}|{alt}" if index % 2 == 0 else f"{alt}|{ref}"
            if bad_gt and donor == "donor_a" and index == 0:
                gt = "A/G"
            rows.append(
                {
                    "chrom": "chr1",
                    "pos": index + 1,
                    "ref": ref,
                    "alt": alt,
                    "GT": gt,
                    "capeak_id": f"peak_{index // 2}",
                    "ref_count": 8,
                    "alt_count": 4,
                    "sample": donor,
                    "is_caqtl": index // 2 % 2 == 0,
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


@pytest.mark.unit
def test_per_donor_peak_analysis_normalizes_phase_and_isolates_models(tmp_path, monkeypatch):
    counts = tmp_path / "counts.tsv"
    _write_peak_counts(counts)
    calls = []

    def fake_fit(path, **kwargs):
        frame = pd.read_csv(path, sep="\t")
        donor = ",".join(frame["donor_id"].drop_duplicates())
        return {
            "fit_id": f"linear:{donor}",
            "method": "linear",
            "rho": None,
            "linear_d1": -3.0,
            "linear_d2": -0.04,
            "n_observations": len(frame),
        }

    def fake_rust(path, *, dispersion_fit, **kwargs):
        frame = pd.read_csv(path, sep="\t")
        assert set(frame["GT"]) == {"0|1", "1|0"}
        assert kwargs["region_col"] == "capeak_id"
        assert kwargs["method"] == "linear"
        calls.append((len(frame), kwargs["phased"]))
        results = []
        for region, group in frame.groupby("capeak_id", sort=True):
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
                    "pval": 0.02,
                    "fdr_pval": 0.04,
                }
            )
        return {
            "results": results,
            "n_observations": len(frame),
            "requested_phased": True,
            "effective_phased": True,
            "dispersion_fit": dispersion_fit,
        }

    monkeypatch.setattr(peak_analysis, "rust_fit_imbalance_dispersion", fake_fit)
    monkeypatch.setattr(peak_analysis, "rust_analyze_imbalance_run", fake_rust)
    paths = peak_analysis.run_per_donor_analysis(
        counts,
        tmp_path / "peak_linear_phased.tsv",
        unit="feature",
        region_col="capeak_id",
        model="linear",
        phased=True,
    )

    assert calls == [(50, True), (50, True)]
    results = pd.read_csv(paths["results"], sep="\t")
    assert results.groupby("donor_id").size().to_dict() == {"donor_a": 25, "donor_b": 25}
    assert results["significant_q05"].all()
    dispersion = pd.read_csv(paths["dispersion"], sep="\t")
    assert dispersion["rho"].isna().all()
    assert dispersion["linear_d1"].eq(-3.0).all()
    qc = pd.read_csv(paths["qc"], sep="\t").set_index("donor_id")
    assert qc.loc["donor_low", "status"] == "excluded_min_observations"
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["analysis"]["scope"] == "per-donor"
    assert provenance["analysis"]["phased"] is True


@pytest.mark.unit
def test_per_donor_peak_analysis_rejects_nonphased_gt(tmp_path, monkeypatch):
    counts = tmp_path / "counts.tsv"
    _write_peak_counts(counts, bad_gt=True)
    monkeypatch.setattr(peak_analysis, "rust_fit_imbalance_dispersion", lambda *a, **k: None)
    monkeypatch.setattr(peak_analysis, "rust_analyze_imbalance_run", lambda *a, **k: None)
    with pytest.raises(ValueError, match="Phased GT must be exact"):
        peak_analysis.run_per_donor_analysis(
            counts,
            tmp_path / "results.tsv",
            unit="feature",
            region_col="capeak_id",
            phased=True,
        )


@pytest.mark.unit
def test_per_donor_peak_analysis_requires_explicit_region(tmp_path, monkeypatch):
    counts = tmp_path / "counts.tsv"
    _write_peak_counts(counts)
    monkeypatch.setattr(peak_analysis, "rust_fit_imbalance_dispersion", lambda *a, **k: None)
    monkeypatch.setattr(peak_analysis, "rust_analyze_imbalance_run", lambda *a, **k: None)
    with pytest.raises(ValueError, match="requires --region-col"):
        peak_analysis.run_per_donor_analysis(
            counts,
            tmp_path / "results.tsv",
            unit="feature",
            region_col="",
        )
