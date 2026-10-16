"""Production-code parity with outputs frozen from unmodified Python 41354fc.

Frozen fixtures are independent of Rust. Live reference tests also require the
pinned archived git objects.
These tests cover statistical compatibility on synthetic counts, not calibration.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scripts.compare_archive_analysis import (
    bulk_counts,
    compare_frames,
    load_archive,
    single_cell_counts,
)
from typer.testing import CliRunner

from analysis.as_analysis_sc import get_imbalance_sc
from analysis.compare_ai import get_compared_imbalance, get_compared_imbalance_diff_snps
from analysis.run_analysis import run_ai_analysis
from analysis.run_analysis_sc import process_adata_inputs

REFERENCE = json.loads(
    (Path(__file__).parent / "fixtures/archive_python_reference.json").read_text()
)


@pytest.mark.parametrize("scenario", ["base", "overlap", "zeros", "unequal_depth"])
@pytest.mark.parametrize("mode", ["single", "shared", "different"])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
@pytest.mark.parametrize("min_count", [0, 10, 20])
def test_live_single_cell_archive_matrix(scenario, mode, phased, pc, min_count):
    from scipy.sparse import csr_matrix

    historical = single_cell_counts()
    if scenario == "overlap":
        historical.uns["feature"] = pd.concat(
            [historical.uns["feature"], pd.DataFrame({"index": [0, 1], "region": ["overlap"] * 2})],
            ignore_index=True,
        )
    elif scenario in ["zeros", "unequal_depth"]:
        for layer in ["ref", "alt"]:
            values = historical.layers[layer].toarray()
            if scenario == "zeros":
                values[0, :] = 0
            else:
                values[:, 2:] *= 4
            historical.layers[layer] = csr_matrix(values)
    # Explicit reference repair: current sums are uint64 to prevent overflow.
    # The unchanged archive remains covered by the frozen tests below.
    archive = load_archive(widen_single_cell_sums=True)
    reference = {
        "single": archive["as_analysis_sc"].get_imbalance_sc,
        "shared": archive["compare_ai"].get_compared_imbalance,
        "different": archive["compare_ai"].get_compared_imbalance_diff_snps,
    }[mode]
    current = {
        "single": get_imbalance_sc,
        "shared": get_compared_imbalance,
        "different": get_compared_imbalance_diff_snps,
    }[mode]
    options = {
        "sample": "sample1",
        "groups": ["A", "B"],
        "phased": phased,
        "pseudocount": pc,
        "min_count": min_count,
    }
    expected = reference(historical.copy(), **options)
    actual = current(historical.T.copy(), **options)
    assert expected.keys() == actual.keys()
    for key in expected:
        comparison = compare_frames(expected[key], actual[key])
        assert comparison["status"] == "PASS", comparison


@pytest.mark.parametrize("phased", [False, True, None])
def test_automatic_sample_selection_filters_like_explicit_selection(phased):
    current = single_cell_counts().T.copy()
    current.var.loc[current.var.index[0], "sample1"] = "0/0"
    automatic = process_adata_inputs(current.copy(), phased=phased)
    explicit = process_adata_inputs(current.copy(), sample="sample1", phased=phased)
    assert automatic.sample == explicit.sample == "sample1"
    assert automatic.phased == explicit.phased
    assert automatic.adata.n_vars == explicit.adata.n_vars == 5
    pd.testing.assert_frame_equal(automatic.adata.var, explicit.adata.var)
    pd.testing.assert_frame_equal(automatic.adata.uns["feature"], explicit.adata.uns["feature"])
    np.testing.assert_array_equal(
        automatic.adata.layers["ref"].toarray(), explicit.adata.layers["ref"].toarray()
    )


@pytest.mark.parametrize("sample", [None, "sample1"])
@pytest.mark.parametrize("requested_phase", [None, True, False])
@pytest.mark.parametrize(
    "genotypes,expected_phase",
    [
        (["0|1", "1|0"] * 3, True),
        (["0|1", "1|0", "0/1", "1/0", "0|1", "1|0"], False),
        (["0|1"] * 6, False),
        (["0/1", "1/0"] * 3, False),
    ],
)
def test_single_cell_phasing_requires_exact_archived_genotype_set(
    sample, requested_phase, genotypes, expected_phase
):
    current = single_cell_counts().T.copy()
    current.var["sample1"] = genotypes
    if requested_phase is not False and not expected_phase:
        with pytest.warns(UserWarning, match="Switching to unphased model"):
            prepared = process_adata_inputs(current, sample=sample, phased=requested_phase)
    else:
        prepared = process_adata_inputs(current, sample=sample, phased=requested_phase)
    assert prepared.phased is (expected_phase and requested_phase is not False)
    assert prepared.adata.n_vars == len(genotypes)


@pytest.mark.parametrize("sample", [None, "sample1"])
@pytest.mark.parametrize("mode", ["single", "shared", "different"])
def test_mixed_phasing_matches_unphased_archived_inference(sample, mode):
    historical = single_cell_counts()
    historical.obs["sample1"] = ["0|1", "1|0", "0/1", "1/0", "0|1", "1|0"]
    with pytest.warns(UserWarning, match="Switching to unphased model"):
        prepared = process_adata_inputs(historical.T.copy(), sample=sample, phased=True)
    archive = load_archive(widen_single_cell_sums=True)
    reference = {
        "single": archive["as_analysis_sc"].get_imbalance_sc,
        "shared": archive["compare_ai"].get_compared_imbalance,
        "different": archive["compare_ai"].get_compared_imbalance_diff_snps,
    }[mode]
    function = {
        "single": get_imbalance_sc,
        "shared": get_compared_imbalance,
        "different": get_compared_imbalance_diff_snps,
    }[mode]
    options = {"sample": "sample1", "groups": prepared.groups, "min_count": 10}
    expected = reference(historical, phased=False, **options)
    actual = function(prepared.adata, phased=prepared.phased, **options)
    assert expected.keys() == actual.keys()
    for key in expected:
        comparison = compare_frames(expected[key], actual[key])
        assert comparison["status"] == "PASS", comparison


@pytest.mark.parametrize("sample", [None, "sample1"])
def test_sample_selection_rejects_no_heterozygotes(sample):
    current = single_cell_counts().T.copy()
    current.var["sample1"] = "0/0"
    with pytest.raises(ValueError, match="Heterozygous genotypes not found"):
        process_adata_inputs(current, sample=sample)


@pytest.mark.parametrize("cutoff", [0.5, 1, 3])
def test_outlier_filtering_matches_archived_feature_reindexing(cutoff):
    from analysis.as_analysis_sc import adata_count_qc

    historical = single_cell_counts()
    historical.obs.loc[historical.obs.index[0], "ref_count"] = 10000
    expected = load_archive()["as_analysis_sc"].adata_count_qc(historical.copy(), z_cutoff=cutoff)
    actual = adata_count_qc(historical.T.copy(), z_cutoff=cutoff)
    np.testing.assert_array_equal(
        actual.layers["ref"].toarray(), expected.layers["ref"].toarray().T
    )
    pd.testing.assert_frame_equal(actual.uns["feature"], expected.uns["feature"])


@pytest.mark.parametrize("phased", [False, True, None])
def test_missing_genotypes_force_unphased_and_missing_groups_fail(phased):
    current = single_cell_counts().T.copy()
    del current.uns["samples"]
    prepared = process_adata_inputs(current.copy(), phased=phased)
    assert prepared.phased is False
    assert prepared.sample is None
    with pytest.raises(KeyError, match="not found"):
        process_adata_inputs(current, groups=["missing"], phased=phased)


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


@pytest.mark.parametrize("command", ["find-imbalance-sc", "compare-imbalance"])
@pytest.mark.parametrize("explicit_sample", [False, True])
def test_single_cell_cli_mixed_phasing_selects_unphased_inference(
    tmp_path, command, explicit_sample
):
    from analysis.__main__ import app

    current = single_cell_counts().T.copy()
    current.var["sample1"] = ["0|1", "1|0", "0/1", "1/0", "0|1", "1|0"]
    counts, groups = tmp_path / "counts.h5ad", tmp_path / "groups.tsv"
    current.write_h5ad(counts)
    groups.write_text("cell0\tA\ncell1\tA\ncell2\tB\ncell3\tB\n")
    arguments = [command, str(counts), str(groups), "--phased", "-o", str(tmp_path / "ai.tsv")]
    if explicit_sample:
        arguments.extend(["--sample", "sample1"])
    with pytest.warns(UserWarning, match="Switching to unphased model"):
        result = CliRunner().invoke(app, arguments)
    assert result.exit_code == 0, result.output + repr(result.exception)
    files = sorted(tmp_path.glob("ai_*.tsv"))
    assert len(files) == (2 if command == "find-imbalance-sc" else 1)
    assert all(not pd.read_csv(file, sep="\t").empty for file in files)
