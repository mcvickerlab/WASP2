import json
from pathlib import Path

import pandas as pd
import pytest

from analysis import run_analysis_per_donor as cohort_analysis
from analysis.run_analysis import run_ai_analysis


def _write_counts(path: Path) -> None:
    rows = []
    for sample, count in [("donor_a", 50), ("donor_b", 51), ("donor_low", 2)]:
        for index in range(count):
            rows.append(
                {
                    "chrom": "chr1",
                    "pos": index + 1,
                    "ref": "A",
                    "alt": "G",
                    "ref_count": 8,
                    "alt_count": 4,
                    "sample": sample,
                    "snv_id": f"chr1:{index + 1}:A:G",
                    "capeak_id": f"peak_{index // 2}",
                    "is_caqtl": index % 2 == 0,
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


@pytest.mark.unit
def test_per_donor_analysis_isolates_rho_and_bh(tmp_path, monkeypatch):
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    calls = []

    def fake_fit(path, **kwargs):
        frame = pd.read_csv(path, sep="\t")
        donor = ",".join(frame["donor_id"].drop_duplicates())
        return {
            "fit_id": f"single:{donor}",
            "method": "single",
            "rho": len(frame) / 1000,
            "linear_d1": None,
            "linear_d2": None,
            "n_observations": len(frame),
        }

    def fake_rust(path, *, dispersion_fit, **kwargs):
        frame = pd.read_csv(path, sep="\t")
        assert "sample" not in frame.columns
        assert "donor_id" not in frame.columns
        calls.append(len(frame))
        results = []
        for row in frame.itertuples(index=False):
            results.append(
                {
                    "region": f"{row.chrom}_{row.pos}",
                    "ref_count": row.ref_count,
                    "alt_count": row.alt_count,
                    "N": row.ref_count + row.alt_count,
                    "snp_count": 1,
                    "null_ll": -2.0,
                    "alt_ll": -1.0,
                    "mu": 0.6,
                    "lrt": 2.0,
                    "pval": 0.02,
                    "fdr_pval": 0.04,
                }
            )
        return {
            "results": results,
            "n_observations": len(frame),
            "requested_phased": False,
            "effective_phased": False,
            "dispersion_fit": dispersion_fit,
        }

    monkeypatch.setattr(cohort_analysis, "rust_fit_imbalance_dispersion", fake_fit)
    monkeypatch.setattr(cohort_analysis, "rust_analyze_imbalance_run", fake_rust)
    output = tmp_path / "results.tsv"
    paths = cohort_analysis.run_per_donor_analysis(
        counts,
        output,
        unit="snv",
        min_donor_observations=50,
    )

    assert sorted(calls) == [50, 51]
    results = pd.read_csv(paths["results"], sep="\t")
    assert results.groupby("donor_id").size().to_dict() == {"donor_a": 50, "donor_b": 51}
    assert results["significant_q05"].all()
    qc = pd.read_csv(paths["qc"], sep="\t").set_index("donor_id")
    assert qc.loc["donor_low", "status"] == "excluded_min_observations"
    assert not bool(qc.loc["donor_low", "analyzed"])
    dispersion = pd.read_csv(paths["dispersion"], sep="\t")
    assert dispersion["donor_id"].tolist() == ["donor_a", "donor_b"]
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["analysis"]["scope"] == "per-donor"
    assert provenance["analysis"]["multiple_testing"].endswith("within each donor")
    assert provenance["donors"] == {"excluded": 1, "included": 2, "total": 3}


@pytest.mark.unit
def test_per_donor_analysis_rejects_conflicting_site_rows(tmp_path, monkeypatch):
    counts = tmp_path / "counts.tsv"
    pd.DataFrame(
        [
            ["chr1", 1, "A", "G", 8, 4, "donor_a"],
            ["chr1", 1, "A", "T", 8, 4, "donor_a"],
        ],
        columns=["chrom", "pos", "ref", "alt", "ref_count", "alt_count", "sample"],
    ).to_csv(counts, sep="\t", index=False)
    monkeypatch.setattr(cohort_analysis, "rust_fit_imbalance_dispersion", lambda *a, **k: None)
    monkeypatch.setattr(cohort_analysis, "rust_analyze_imbalance_run", lambda *a, **k: None)
    with pytest.raises(ValueError, match="Conflicting alleles"):
        cohort_analysis.run_per_donor_analysis(
            counts,
            tmp_path / "results.tsv",
            unit="snv",
            min_donor_observations=1,
        )


@pytest.mark.unit
def test_per_donor_analysis_enforces_manifest_hash_and_contract(tmp_path, monkeypatch):
    counts = tmp_path / "counts.tsv"
    _write_counts(counts)
    monkeypatch.setattr(cohort_analysis, "rust_fit_imbalance_dispersion", lambda *a, **k: None)
    monkeypatch.setattr(cohort_analysis, "rust_analyze_imbalance_run", lambda *a, **k: None)
    with pytest.raises(ValueError, match="requires a locked count bundle"):
        cohort_analysis.run_per_donor_analysis(
            counts,
            tmp_path / "results.tsv",
            unit="snv",
            expected_manifest_sha256="0" * 64,
        )
    with pytest.raises(ValueError, match="unphased"):
        run_ai_analysis(counts, scope="per-donor", unit="snv", phased=True)
