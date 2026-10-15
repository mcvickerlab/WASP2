"""Repeat the full 32-run real ATAC matrix against executed archived Python.

Count tables are frozen measurements from public reads, not Rust-generated
expectations. Linux Python 3.11 CI pins the original numerical environment.
"""

import contextlib
import gzip
import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy
from scipy.special import expit
from scripts.compare_archive_analysis import adapt_linear_reference, compare_frames, load_archive

from analysis.run_analysis_per_donor import run_per_donor_analysis

FIXTURES = Path(__file__).parent / "fixtures/public_atac"
METADATA = json.loads((FIXTURES / "metadata.json").read_text())


@pytest.mark.rust
@pytest.mark.skipif(
    (np.__version__, scipy.__version__) != ("1.26.4", "1.13.1"),
    reason="Full historical matrix uses the pinned NumPy/SciPy environment (Python 3.11 CI)",
)
@pytest.mark.parametrize("unit", ["snv", "peak"])
@pytest.mark.parametrize("model", ["single", "linear"])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
def test_full_public_atac_matrix(tmp_path, unit, model, phased, pc):
    import wasp2_rust

    source = METADATA["counts"][unit]
    raw = gzip.decompress((FIXTURES / source["file"]).read_bytes())
    assert hashlib.sha256(raw).hexdigest() == source["sha256"]
    frame = pd.read_csv(io.BytesIO(raw), sep="\t")
    assert len(frame) == source["rows"]
    if phased:
        frame = frame.loc[frame["GT"].isin(["0|1", "1|0"])].copy()
    region = "region" if unit == "peak" else None
    modules = load_archive(METADATA["archive_commit"])
    if model == "linear":
        adapt_linear_reference(modules)
    archive = modules["as_analysis"]
    fitted = []
    optimizer_name = "minimize_scalar" if model == "single" else "minimize"
    original_optimizer = getattr(archive, optimizer_name)

    def capture_fit(*args, **kwargs):
        result = original_optimizer(*args, **kwargs)
        if not fitted:
            fitted.append(result.x)
        return result

    setattr(archive, optimizer_name, capture_fit)
    if model == "linear" and phased:
        # Explicit reference repair: archived linear_model ignores this flag.
        original_parse = archive.parse_opt
        archive.parse_opt = lambda df, disp=None, phased=False: original_parse(
            df, disp=disp, phased=True
        )
    with contextlib.redirect_stdout(io.StringIO()):
        expected = archive.get_imbalance(
            frame.copy(), method=model, phased=phased, region_col=region, pseudocount=pc
        )
    expected = expected.rename(columns={"variant": "region"})
    path = tmp_path / "counts.tsv"
    frame.to_csv(path, sep="\t", index=False)
    output = wasp2_rust.analyze_imbalance_run(
        str(path), method=model, phased=phased, region_col=region, pseudocount=pc
    )
    comparison = compare_frames(expected, pd.DataFrame(output["results"]))
    assert comparison["status"] == "PASS", json.dumps(
        {
            "comparison": comparison,
            "python_fit": np.asarray(fitted[0]).tolist(),
            "rust_fit": {key: output[key] for key in ["rho", "linear_d1", "linear_d2"]},
        }
    )

    def check_dispersion(fit):
        if model == "single":
            np.testing.assert_allclose(fit["rho"], fitted[0], rtol=5e-5, atol=5e-5)
        else:
            depth = (frame["ref_count"] + frame["alt_count"]).to_numpy()
            depth = depth[depth >= 10]
            actual = expit(
                fit["linear_d1"]
                + ((depth - fit["linear_depth_center"]) / fit["linear_depth_scale"])
                * fit["linear_d2"]
            )
            expected_rho = expit(fitted[0][0] + (depth + 2 * pc) * fitted[0][1])
            np.testing.assert_allclose(actual, expected_rho, rtol=5e-5, atol=5e-5)

    check_dispersion(output)
    cli_path = tmp_path / "cli.tsv"
    command = [
        sys.executable,
        "-c",
        "from analysis.__main__ import app; app()",
        "find-imbalance",
        str(path),
        "--model",
        model,
        "--pseudocount",
        str(pc),
        "--out_file",
        str(cli_path),
    ]
    command += ["--region_col", "region"] if region else ["--per-variant"]
    if phased:
        command.append("--phased")
    subprocess.run(command, check=True, capture_output=True, text=True)
    comparison = compare_frames(expected, pd.read_csv(cli_path, sep="\t"))
    assert comparison["status"] == "PASS", ("CLI", comparison)
    if unit == "peak":
        frame.assign(donor_id="HG001").to_csv(path, sep="\t", index=False)
        for scope in ["global", "per-donor"]:
            paths = run_per_donor_analysis(
                path,
                tmp_path / f"{scope}.tsv",
                unit="feature",
                model=model,
                region_col="region",
                phased=phased,
                pseudocount=pc,
                dispersion_scope=scope,
            )
            actual = pd.read_csv(paths["results"], sep="\t").rename(
                columns={"feature_id": "region"}
            )
            comparison = compare_frames(expected, actual)
            assert comparison["status"] == "PASS", (scope, comparison)
            check_dispersion(pd.read_csv(paths["dispersion"], sep="\t").iloc[0])
