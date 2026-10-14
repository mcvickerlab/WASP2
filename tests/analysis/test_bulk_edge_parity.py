"""Execute archived Python on mixed GT and interior dispersion boundaries."""

import contextlib
import io
import json

import numpy as np
import pandas as pd
import pytest
from scipy.stats import betabinom
from scripts.compare_archive_analysis import bulk_counts, compare_frames, load_archive

from analysis.run_analysis_per_donor import run_per_donor_analysis


@pytest.mark.rust
@pytest.mark.parametrize("gt", ["0|0", "1|1"])
@pytest.mark.parametrize("pc", [0, 1])
@pytest.mark.parametrize("scope", [None, "per-donor", "global"])
def test_mixed_numeric_genotypes_match_python(tmp_path, gt, pc, scope):
    import wasp2_rust

    frame = bulk_counts().drop(columns="gene")
    frame["GT"] = "0|1" if gt == "0|0" else "1|0"
    frame.loc[frame.index[1::2], "GT"] = gt
    with contextlib.redirect_stdout(io.StringIO()):
        expected = load_archive()["as_analysis"].get_imbalance(
            frame.copy(), phased=True, region_col="region", pseudocount=pc
        )
    path = tmp_path / "counts.tsv"
    if scope is None:
        frame.to_csv(path, sep="\t", index=False)
        output = wasp2_rust.analyze_imbalance_run(
            str(path), phased=True, region_col="region", pseudocount=pc
        )
        assert output["effective_phased"] is True
        actual = pd.DataFrame(output["results"])
    else:
        frame.assign(donor_id="A").to_csv(path, sep="\t", index=False)
        paths = run_per_donor_analysis(
            path,
            tmp_path / "result.tsv",
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


@pytest.mark.rust
@pytest.mark.parametrize("gt", ["0/1", "A|G", ".|."])
@pytest.mark.parametrize("scope", [None, "per-donor", "global"])
def test_malformed_mixed_genotypes_fail_like_python(tmp_path, gt, scope):
    import wasp2_rust

    frame = bulk_counts().drop(columns="gene")
    frame.loc[0, "GT"] = gt
    with contextlib.redirect_stdout(io.StringIO()), pytest.raises((ValueError, TypeError)):
        load_archive()["as_analysis"].get_imbalance(frame.copy(), phased=True, region_col="region")
    path = tmp_path / "counts.tsv"
    if scope is None:
        frame.to_csv(path, sep="\t", index=False)
        with pytest.raises(RuntimeError, match="Invalid GT"):
            wasp2_rust.analyze_imbalance_run(str(path), phased=True, region_col="region")
    else:
        frame.assign(donor_id="A").to_csv(path, sep="\t", index=False)
        with pytest.raises((ValueError, TypeError)):
            run_per_donor_analysis(
                path,
                tmp_path / "result.tsv",
                unit="feature",
                region_col="region",
                phased=True,
                dispersion_scope=scope,
                min_donor_observations=4,
            )


@pytest.mark.rust
@pytest.mark.parametrize("scope", ["per-donor", "global"])
@pytest.mark.parametrize("orientation", ["0|1", "1|0"])
@pytest.mark.parametrize("model", ["single", "linear"])
@pytest.mark.parametrize("pc", [0, 1])
def test_numeric_gt_duplicates_preserve_pre_dedup_phase_decision(
    tmp_path, scope, orientation, model, pc
):
    frame = bulk_counts().drop(columns="gene")
    frame["GT"] = orientation
    duplicate = frame.iloc[:1].copy()
    duplicate["GT"] = "0|0" if orientation == "0|1" else "1|1"
    frame = pd.concat([frame, duplicate], ignore_index=True)
    from scripts.compare_archive_analysis import adapt_linear_reference

    modules = load_archive()
    if model == "linear":
        adapt_linear_reference(modules)
        original = modules["as_analysis"].parse_opt
        modules["as_analysis"].parse_opt = lambda df, disp=None, phased=False: original(
            df, disp=disp, phased=True
        )
    with contextlib.redirect_stdout(io.StringIO()):
        expected = modules["as_analysis"].get_imbalance(
            frame.copy(),
            phased=True,
            region_col="region",
            method=model,
            pseudocount=pc,
        )
    path = tmp_path / "counts.tsv"
    frame.assign(donor_id="A").to_csv(path, sep="\t", index=False)
    paths = run_per_donor_analysis(
        path,
        tmp_path / "result.tsv",
        unit="feature",
        region_col="region",
        phased=True,
        model=model,
        pseudocount=pc,
        dispersion_scope=scope,
        min_donor_observations=4,
    )
    actual = pd.read_csv(paths["results"], sep="\t").rename(columns={"feature_id": "region"})
    comparison = compare_frames(expected, actual)
    assert comparison["status"] == "PASS", comparison
    fit = pd.read_csv(paths["dispersion"], sep="\t").iloc[0]
    assert fit["fit_observations"] == len(frame) - 1
    qc = pd.read_csv(paths["qc"], sep="\t").iloc[0]
    assert qc["input_rows"] == len(frame)
    assert qc["eligible_test_rows"] == len(frame) - 1
    provenance = json.loads(paths["provenance"].read_text())
    assert provenance["donors"]["phasing"][0]["effective_phased"] is True


@pytest.mark.rust
@pytest.mark.parametrize("rho", [1e-12, 1e-10, 1e-8, 1.0 - 1e-12])
@pytest.mark.parametrize("ref_count", [0, 10, 20])
def test_interior_dispersion_is_not_clipped(tmp_path, rho, ref_count):
    import wasp2_rust

    path = tmp_path / "counts.tsv"
    pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "pos": 101,
                "ref": "A",
                "alt": "G",
                "region": "peak",
                "ref_count": ref_count,
                "alt_count": 20 - ref_count,
            }
        ]
    ).to_csv(path, sep="\t", index=False)
    output = wasp2_rust.analyze_imbalance_run(
        str(path), pseudocount=0, region_col="region", rho=rho
    )
    shape = 0.5 * (1.0 - rho) / rho
    expected = betabinom.logpmf(ref_count, 20, shape, shape)
    assert output["rho"] == rho
    np.testing.assert_allclose(output["results"][0]["null_ll"], expected, rtol=5e-5, atol=5e-5)
