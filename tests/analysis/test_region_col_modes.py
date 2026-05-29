"""Regression tests for region_col grouping modes + CLI flag forwarding.

Covers the two grouping modes selectable via the Rust ``analyze_imbalance(region_col=...)``
flag and the ``wasp2 analysis`` CLI plumbing in ``analysis.run_analysis``:

* feature/peak mode (``region_col="region"``) groups variants within a region;
* SNV-solo / per-variant mode (``region_col=None``) tests each SNV independently and
  labels regions ``{chrom}_{pos}`` (must match Python's key for parity);
* ``--per-variant`` forces per-variant even when a region column is present;
* ``--per-variant`` + ``--region_col`` is rejected;
* ``--groupby`` fails closed on the Rust backend (it only re-keys the grouping column,
  so ``--region_col <parent>`` is the supported equivalent);
* ``--phased`` / ``--region_col`` are forwarded to the Rust backend.
"""

from __future__ import annotations

import pytest

# --------------------------------------------------------------------------- helpers


def _write_counts(path, *, with_region: bool = True):
    """Write a tiny phased count TSV (two het SNVs in one peak).

    ``with_region=True`` emits the 9-column layout that carries a ``region`` column;
    ``with_region=False`` omits it. Only the header is consulted by ``WaspAnalysisData``.
    """
    if with_region:
        header = "chrom\tstart\tpos\tref\tregion\tGT\tref_count\talt_count\tN\n"
        rows = "chr1\t100\t101\tA\tpeakX\t0|1\t10\t8\t18\nchr1\t200\t201\tC\tpeakX\t1|0\t7\t9\t16\n"
    else:
        header = "chrom\tstart\tpos\tref\tGT\tref_count\talt_count\tN\n"
        rows = "chr1\t100\t101\tA\t0|1\t10\t8\t18\nchr1\t200\t201\tC\t1|0\t7\t9\t16\n"
    path.write_text(header + rows)
    return path


# --------------------------------------------------- Rust backend grouping modes


@pytest.mark.rust
def test_rust_feature_mode_groups_by_region(tmp_path):
    """region_col='region' collapses both SNVs into one peak-level result."""
    wasp2_rust = pytest.importorskip("wasp2_rust")
    tsv = _write_counts(tmp_path / "counts.tsv")
    res = wasp2_rust.analyze_imbalance(
        str(tsv),
        min_count=0,
        pseudocount=0,
        method="single",
        phased=True,
        region_col="region",
    )
    assert len(res) == 1
    assert res[0]["region"] == "peakX"
    assert res[0]["snp_count"] == 2
    assert "mu" in res[0]  # phased model engaged


@pytest.mark.rust
def test_rust_solo_mode_per_variant_labels(tmp_path):
    """region_col=None tests each SNV independently, labelled {chrom}_{pos}."""
    wasp2_rust = pytest.importorskip("wasp2_rust")
    tsv = _write_counts(tmp_path / "counts.tsv")
    res = wasp2_rust.analyze_imbalance(
        str(tsv),
        min_count=0,
        pseudocount=0,
        method="single",
        phased=True,
        region_col=None,
    )
    regions = sorted(r["region"] for r in res)
    # {chrom}_{pos} (matches Python's df["chrom"] + "_" + df["pos"]) — NOT chrom_start_pos.
    assert regions == ["chr1_101", "chr1_201"]
    assert all(r["snp_count"] == 1 for r in res)


@pytest.mark.rust
def test_rust_region_col_missing_raises(tmp_path):
    """A region_col naming a column that is absent must error (fail closed)."""
    wasp2_rust = pytest.importorskip("wasp2_rust")
    tsv = _write_counts(tmp_path / "counts.tsv")
    with pytest.raises(RuntimeError):
        wasp2_rust.analyze_imbalance(
            str(tsv),
            min_count=0,
            pseudocount=0,
            method="single",
            phased=True,
            region_col="does_not_exist",
        )


# --------------------------------------------------- CLI plumbing (no extension)


@pytest.fixture
def spy_rust(monkeypatch):
    """Replace run_analysis.rust_analyze_imbalance with a kwargs-capturing spy.

    Setting it to a non-None callable also bypasses the "extension not available"
    guard, so these tests run without the built Rust extension.
    """
    from analysis import run_analysis

    captured: dict = {}

    def _spy(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        return [{"region": "x", "fdr_pval": 0.5}]  # minimal DataFrame-able result

    monkeypatch.setattr(run_analysis, "rust_analyze_imbalance", _spy)
    return captured


@pytest.mark.unit
def test_cli_forwards_phased_and_region_col(tmp_path, spy_rust):
    from analysis.run_analysis import run_ai_analysis

    tsv = _write_counts(tmp_path / "counts.tsv")
    run_ai_analysis(
        str(tsv),
        min_count=0,
        pseudocount=0,
        phased=True,
        region_col="region",
        out_file=str(tmp_path / "out.tsv"),
    )
    assert spy_rust["kwargs"].get("phased") is True
    assert spy_rust["kwargs"].get("region_col") == "region"


@pytest.mark.unit
def test_cli_groupby_fails_closed(tmp_path, spy_rust):
    from analysis.run_analysis import run_ai_analysis

    tsv = _write_counts(tmp_path / "counts.tsv")
    with pytest.raises(RuntimeError, match="groupby"):
        run_ai_analysis(
            str(tsv),
            min_count=0,
            pseudocount=0,
            groupby="parent",
            out_file=str(tmp_path / "out.tsv"),
        )


@pytest.mark.unit
def test_cli_per_variant_forces_solo(tmp_path, spy_rust):
    from analysis.run_analysis import run_ai_analysis

    # Region column IS present; --per-variant must still force region_col=None.
    tsv = _write_counts(tmp_path / "counts.tsv", with_region=True)
    run_ai_analysis(
        str(tsv),
        min_count=0,
        pseudocount=0,
        phased=True,
        per_variant=True,
        out_file=str(tmp_path / "out.tsv"),
    )
    assert spy_rust["kwargs"].get("region_col") is None


@pytest.mark.unit
def test_per_variant_region_col_conflict(tmp_path):
    from analysis.run_analysis import run_ai_analysis

    tsv = _write_counts(tmp_path / "counts.tsv")
    with pytest.raises(ValueError):
        run_ai_analysis(
            str(tsv),
            min_count=0,
            pseudocount=0,
            per_variant=True,
            region_col="region",
            out_file=str(tmp_path / "out.tsv"),
        )
