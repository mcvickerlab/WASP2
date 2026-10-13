"""Production-code parity with outputs frozen from unmodified Python 41354fc.

The fixture is independent of Rust and does not need archived git objects in CI.
These tests cover statistical compatibility on synthetic counts, not calibration.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scripts.compare_archive_analysis import bulk_counts, compare_frames, single_cell_counts
from typer.testing import CliRunner

from analysis.as_analysis_sc import get_imbalance_sc
from analysis.compare_ai import get_compared_imbalance, get_compared_imbalance_diff_snps
from analysis.run_analysis import run_ai_analysis
from analysis.run_analysis_sc import process_adata_inputs

REFERENCE = json.loads(
    (Path(__file__).parent / "fixtures/archive_python_reference.json").read_text()
)


@pytest.mark.rust
@pytest.mark.parametrize("unit", ["snv", "atac", "rna"])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
def test_bulk_archive_parity(tmp_path, unit, phased, pc):
    region = {"snv": None, "atac": "region", "rna": "gene"}[unit]
    counts = bulk_counts().drop(
        columns=["region", "gene"] if unit == "snv" else ["gene"] if unit == "atac" else ["region"]
    )
    count_file, output = tmp_path / "counts.tsv", tmp_path / "results.tsv"
    counts.to_csv(count_file, sep="\t", index=False)
    run_ai_analysis(
        count_file,
        model="single",
        phased=phased,
        region_col=region,
        pseudocount=pc,
        out_file=str(output),
    )
    expected = pd.DataFrame(REFERENCE["bulk"][f"{unit}/{phased}/{pc}"])
    actual = pd.read_csv(output, sep="\t")
    comparison = compare_frames(expected, actual)
    assert comparison["status"] == "PASS", comparison


@pytest.mark.rust
@pytest.mark.parametrize(
    "unit,phased", [("snv", False), ("snv", True), ("atac", False), ("rna", False)]
)
@pytest.mark.parametrize("pc", [0, 1])
def test_linear_shape_adapted_reference(tmp_path, unit, phased, pc):
    # The fixture records the reference-only adaptations. Phased multi-SNP
    # linear results are excluded because the archive silently ignored phasing.
    region = {"snv": None, "atac": "region", "rna": "gene"}[unit]
    counts = bulk_counts().drop(
        columns=["region", "gene"] if unit == "snv" else ["gene"] if unit == "atac" else ["region"]
    )
    count_file, output = tmp_path / "counts.tsv", tmp_path / "results.tsv"
    counts.to_csv(count_file, sep="\t", index=False)
    run_ai_analysis(
        count_file,
        model="linear",
        phased=phased,
        region_col=region,
        pseudocount=pc,
        out_file=str(output),
    )
    expected = pd.DataFrame(REFERENCE["linear"][f"{unit}/{phased}/{pc}"])
    comparison = compare_frames(expected, pd.read_csv(output, sep="\t"))
    assert comparison["status"] == "PASS", comparison


@pytest.mark.parametrize("mode", ["single", "shared", "different"])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
def test_single_cell_archive_parity(mode, phased, pc):
    current = single_cell_counts().T.copy()
    prepared = process_adata_inputs(current, sample="sample1", phased=phased)
    function = {
        "single": get_imbalance_sc,
        "shared": get_compared_imbalance,
        "different": get_compared_imbalance_diff_snps,
    }[mode]
    actual = function(
        prepared.adata,
        sample=prepared.sample,
        phased=prepared.phased,
        groups=prepared.groups,
        pseudocount=pc,
        min_count=10,
    )
    expected = REFERENCE["single_cell"][f"{mode}/{phased}/{pc}"]
    assert set(expected) == {str(k) for k in actual}
    for key, frame in actual.items():
        comparison = compare_frames(pd.DataFrame(expected[str(key)]), frame)
        assert comparison["status"] == "PASS", comparison


def test_barcode_mapping_and_genotype_filter_preserve_feature_indices(tmp_path):
    current = single_cell_counts().T.copy()
    current.var.loc[current.var.index[1], "sample1"] = "0/0"
    mapping = tmp_path / "groups.tsv"
    mapping.write_text("cell0\tnewA\ncell1\tnewA\ncell2\tnewB\ncell3\tnewB\n")
    prepared = process_adata_inputs(current, bc_map=mapping, sample="sample1", phased=False)
    assert prepared.adata.shape == (4, 5)
    assert prepared.adata.obs["group"].tolist() == ["newA", "newA", "newB", "newB"]
    assert prepared.adata.var["index"].tolist() == list(range(5))
    assert prepared.adata.uns["feature"]["index"].tolist() == list(range(5))
    assert prepared.adata.uns["feature"]["region"].tolist() == ["p0", "p1", "p2", "p2", "p2"]


def test_unphased_counts_without_genotype_metadata():
    current = single_cell_counts().T.copy()
    del current.uns["samples"]
    current.var = current.var.drop(columns="sample1")
    prepared = process_adata_inputs(current)
    assert prepared.sample is None
    assert prepared.phased is False
    assert set(get_imbalance_sc(prepared.adata, groups=prepared.groups)) == {"A", "B"}


@pytest.mark.parametrize("mode", ["shared", "different"])
def test_differential_sums_do_not_wrap_at_uint16(monkeypatch, mode):
    import analysis.compare_ai as comparison

    current = single_cell_counts().T.copy()
    current.layers["ref"] = current.layers["ref"].astype(np.uint32) * 10000
    current.layers["alt"] = current.layers["alt"].astype(np.uint32) * 10000
    captured = []
    if mode == "shared":

        def capture(disp, ref1, n1, ref2, n2, regions, gt):
            captured.append((ref1.copy(), n1.copy(), ref2.copy(), n2.copy()))
            return pd.DataFrame({"region": ["p0"]})

        monkeypatch.setattr(comparison, "compare_imbalance_between_groups", capture)
        function = comparison.get_compared_imbalance
    else:

        def capture(disp, ref1, n1, gt1, regions1, ref2, n2, gt2, regions2):
            captured.append((ref1.copy(), n1.copy(), ref2.copy(), n2.copy()))
            return pd.DataFrame({"region": ["p0"]})

        monkeypatch.setattr(comparison, "compare_imbalance_between_groups_diff_snps", capture)
        function = comparison.get_compared_imbalance_diff_snps
    function(current, sample="sample1", phased=False, pseudocount=0)
    assert len(captured) == 1
    ref1, n1, ref2, n2 = captured[0]
    np.testing.assert_array_equal(ref1, np.asarray(current.layers["ref"][:2].sum(axis=0)).ravel())
    np.testing.assert_array_equal(ref2, np.asarray(current.layers["ref"][2:].sum(axis=0)).ravel())
    np.testing.assert_array_equal(n1, np.full(6, 200000))
    np.testing.assert_array_equal(n2, np.full(6, 200000))


@pytest.mark.parametrize("command", ["find-imbalance-sc", "compare-imbalance"])
def test_single_cell_cli_reads_current_count_layout(tmp_path, command):
    from analysis.__main__ import app

    counts, groups, output = tmp_path / "counts.h5ad", tmp_path / "groups.tsv", tmp_path / "ai.tsv"
    single_cell_counts().T.copy().write_h5ad(counts)
    groups.write_text("cell0\tA\ncell1\tA\ncell2\tB\ncell3\tB\n")
    result = CliRunner().invoke(
        app, [command, str(counts), str(groups), "--sample", "sample1", "--out_file", str(output)]
    )
    assert result.exit_code == 0, result.output + repr(result.exception)
    files = sorted(tmp_path.glob("ai_*.tsv"))
    assert len(files) == (2 if command == "find-imbalance-sc" else 1)
    assert all(not pd.read_csv(file, sep="\t").empty for file in files)
