#!/usr/bin/env python3
"""Compare analysis with Aaron Ho's pinned, pre-migration Python implementation.

Run from a full repository clone with an installed Rust extension:
    python scripts/compare_archive_analysis.py --output comparison.json
The archived code is read from git without altering it. AnnData is transposed
at the boundary to compare the historical SNP-by-cell layout with today's
cell-by-SNP layout. This is a compatibility check, not statistical calibration.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import subprocess
import sys
import tempfile
import types
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_REF = "41354fc"
STAT_COLUMNS = ["null_ll", "alt_ll", "mu", "lrt", "pval", "fdr_pval"]
COUNT_COLUMNS = ["ref_count", "alt_count", "N", "snp_count", "num_snps"]


def load_archive(ref: str = ARCHIVE_REF) -> dict[str, types.ModuleType]:
    """Load unmodified historical modules, restoring their legacy import names."""
    names = ["as_analysis", "as_analysis_sc", "compare_ai"]
    saved = {name: sys.modules.get(name) for name in names}
    modules = {}
    try:
        for name in names:
            source = subprocess.check_output(
                ["git", "show", f"{ref}:src/analysis/{name}.py"], cwd=ROOT, text=True
            )
            module = types.ModuleType(f"archive_{name}")
            modules[name] = module
            sys.modules[name] = module
            exec(compile(source, f"{ref}:src/analysis/{name}.py", "exec"), module.__dict__)
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
    return modules


def adapt_linear_reference(archive: dict[str, types.ModuleType]) -> None:
    """Explicitly repair archived array shapes, without changing its likelihood.

    The archive passes length-one arrays where SciPy requires scalar objectives,
    and passes the full dispersion vector to the first/subsequent SNP slices.
    This optional reference-only adapter aligns those shapes. It deliberately
    retains the archived optimizer, first-SNP anchor, chi-square p-values, and
    linear model's failure to forward its phased argument.
    """
    module = archive["as_analysis"]
    original_prob = module.opt_prob

    def scalar_prob(*args, **kwargs):
        result = original_prob(*args, **kwargs)
        return float(np.asarray(result).item()) if np.asarray(result).size == 1 else result

    def aligned_dp(prob, disp, first_ref, first_n, phase_ref, phase_n):
        dispersions = np.asarray(disp)
        per_snp = dispersions.ndim and len(dispersions) > 1
        first_disp = dispersions[0] if per_snp else disp
        phase_disp = dispersions[1:] if per_snp else disp
        first_ll = scalar_prob(prob, first_disp, first_ref[0], first_n[0])
        p1 = np.atleast_1d(scalar_prob(prob, phase_disp, phase_ref, phase_n, log=False))
        p2 = np.atleast_1d(scalar_prob(1 - prob, phase_disp, phase_ref, phase_n, log=False))
        return float(first_ll - np.log(np.prod(0.5 * p1 + 0.5 * p2)))

    module.opt_prob = scalar_prob
    module.opt_unphased_dp = aligned_dp


def bulk_counts() -> pd.DataFrame:
    pairs = [
        (18, 2),
        (3, 17),
        (22, 18),
        (25, 15),
        (29, 11),
        (10, 30),
        (10, 10),
        (8, 12),
        (35, 5),
        (6, 34),
        (14, 6),
        (7, 13),
    ]
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * len(pairs),
            "pos": np.arange(1, len(pairs) + 1) * 100,
            "ref": ["A"] * len(pairs),
            "alt": ["G"] * len(pairs),
            "GT": ["0|1", "1|0"] * (len(pairs) // 2),
            "region": [f"peak{i // 2}" for i in range(len(pairs))],
            "gene": [f"gene{i // 3}" for i in range(len(pairs))],
            "ref_count": [p[0] for p in pairs],
            "alt_count": [p[1] for p in pairs],
            "other_count": [0] * len(pairs),
        }
    )


def single_cell_counts() -> AnnData:
    """Historical layout; includes both single- and multi-SNP features."""
    ref = np.array(
        [[9, 8, 2, 3], [2, 3, 8, 9], [6, 5, 5, 4], [8, 7, 7, 8], [3, 2, 2, 3], [5, 4, 6, 5]],
        dtype=np.uint16,
    )
    alt = 10 - ref
    obs = pd.DataFrame(
        {
            "index": np.arange(6),
            "sample1": ["0|1", "1|0"] * 3,
            "ref_count": ref.sum(axis=1),
            "alt_count": alt.sum(axis=1),
        },
        index=[f"snp{i}" for i in range(6)],
    )
    var = pd.DataFrame({"group": ["A", "A", "B", "B"]}, index=[f"cell{i}" for i in range(4)])
    adata = AnnData(X=csr_matrix(ref + alt), obs=obs, var=var)
    adata.layers["ref"] = csr_matrix(ref)
    adata.layers["alt"] = csr_matrix(alt)
    adata.uns["samples"] = ["sample1"]
    adata.uns["feature"] = pd.DataFrame(
        {"index": np.arange(6), "region": ["p0", "p0", "p1", "p2", "p2", "p2"]}
    )
    return adata


def compare_frames(expected: pd.DataFrame, observed: pd.DataFrame, key: str = "region") -> dict:
    expected = expected.sort_values(key).reset_index(drop=True)
    observed = observed.sort_values(key).reset_index(drop=True)
    if expected[key].tolist() != observed[key].tolist():
        return {"status": "FAIL", "reason": "different tested identities"}
    deltas = {}
    failures = []
    for column in [*COUNT_COLUMNS, *STAT_COLUMNS, "mu1", "mu2", "combined_mu"]:
        if column not in expected:
            continue
        if column not in observed:
            failures.append(f"missing {column}")
            continue
        left = expected[column].to_numpy(dtype=float)
        right = observed[column].to_numpy(dtype=float)
        deltas[column] = float(np.max(np.abs(left - right))) if len(left) else 0.0
        if column in COUNT_COLUMNS:
            equal = np.array_equal(left, right)
        else:
            equal = (
                np.all(np.isfinite(left))
                and np.all(np.isfinite(right))
                and np.allclose(left, right, rtol=5e-5, atol=5e-5)
            )
        if not equal:
            failures.append(column)
    for alpha in [0.05, 0.1]:
        if not np.array_equal(expected["fdr_pval"] < alpha, observed["fdr_pval"] < alpha):
            failures.append(f"significance at q={alpha}")
    return {
        "status": "FAIL" if failures else "PASS",
        "failures": failures,
        "max_absolute_differences": deltas,
    }


def run_comparison(archive: dict[str, types.ModuleType]) -> dict:
    import wasp2_rust
    from analysis.as_analysis_sc import get_imbalance_sc
    from analysis.compare_ai import get_compared_imbalance, get_compared_imbalance_diff_snps
    from analysis.run_analysis import run_ai_analysis
    from analysis.run_analysis_sc import process_adata_inputs

    cases = {}
    with tempfile.TemporaryDirectory(prefix="wasp2-parity-") as tmp:
        temp = Path(tmp)
        for model in ["single", "linear"]:
            for unit in ["snv", "atac", "rna"]:
                for phased in [False, True]:
                    for pc in [0, 1]:
                        name = f"bulk/{model}/{unit}/phased={phased}/pc={pc}"
                        region = {"snv": None, "atac": "region", "rna": "gene"}[unit]
                        counts = bulk_counts().drop(
                            columns=["region", "gene"]
                            if unit == "snv"
                            else ["gene"]
                            if unit == "atac"
                            else ["region"]
                        )
                        counts.to_csv(temp / "counts.tsv", sep="\t", index=False)
                        try:
                            expected = archive["as_analysis"].get_imbalance(
                                counts.copy(),
                                method=model,
                                phased=phased,
                                region_col=region,
                                pseudocount=pc,
                            )
                        except Exception as exc:
                            cases[name] = {"status": "REFERENCE_ERROR", "reason": repr(exc)}
                            continue
                        try:
                            run_ai_analysis(
                                temp / "counts.tsv",
                                model=model,
                                phased=phased,
                                region_col=region,
                                pseudocount=pc,
                                out_file=str(temp / "results.tsv"),
                            )
                            actual = pd.read_csv(temp / "results.tsv", sep="\t")
                            if region != "region":
                                expected = expected.rename(columns={region or "variant": "region"})
                            cases[name] = compare_frames(expected, actual)
                        except Exception as exc:
                            cases[name] = {"status": "ERROR", "reason": repr(exc)}
        for phased in [False, True]:
            for pc in [0, 1]:
                old = single_cell_counts()
                new = old.T.copy()
                options = {
                    "sample": "sample1",
                    "phased": phased,
                    "pseudocount": pc,
                    "min_count": 10,
                }
                for label, reference, current in [
                    ("single-cell", archive["as_analysis_sc"].get_imbalance_sc, get_imbalance_sc),
                    (
                        "differential/shared",
                        archive["compare_ai"].get_compared_imbalance,
                        get_compared_imbalance,
                    ),
                    (
                        "differential/different",
                        archive["compare_ai"].get_compared_imbalance_diff_snps,
                        get_compared_imbalance_diff_snps,
                    ),
                ]:
                    name = f"{label}/phased={phased}/pc={pc}"
                    try:
                        expected = reference(old.copy(), **options)
                    except Exception as exc:
                        cases[name] = {"status": "REFERENCE_ERROR", "reason": repr(exc)}
                        continue
                    try:
                        prepared = process_adata_inputs(new.copy(), sample="sample1", phased=phased)
                        actual = current(prepared.adata, **options)
                        details = {
                            str(key): compare_frames(frame, actual[key])
                            for key, frame in expected.items()
                        }
                        cases[name] = {
                            "status": "PASS"
                            if expected.keys() == actual.keys()
                            and all(d["status"] == "PASS" for d in details.values())
                            else "FAIL",
                            "groups": details,
                        }
                    except Exception as exc:
                        cases[name] = {"status": "ERROR", "reason": repr(exc)}
    return {
        "archive_commit": subprocess.check_output(
            ["git", "rev-parse", ARCHIVE_REF], cwd=ROOT, text=True
        ).strip(),
        "current_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "extension": wasp2_rust.__file__,
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--adapt-linear-reference",
        action="store_true",
        help="Explicitly repair archived scalar/array shapes; report this adaptation in JSON.",
    )
    args = parser.parse_args()
    # Import the installed extension before src's stub package can shadow it.
    import wasp2_rust  # noqa: F401

    sys.path.insert(0, str(ROOT / "src"))
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore", FutureWarning)
        archive = load_archive()
        if args.adapt_linear_reference:
            adapt_linear_reference(archive)
        report = run_comparison(archive)
    report["reference_adaptations"] = (
        ["Scalarize one-element objectives and align per-SNP dispersion arrays"]
        if args.adapt_linear_reference
        else []
    )
    report["known_reference_limitations"] = [
        "Archived linear_model does not forward phased to parse_opt; phased multi-SNP "
        "linear comparisons therefore test different declared behavior."
    ]
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    summary = {
        status: sum(c["status"] == status for c in report["cases"].values())
        for status in ["PASS", "FAIL", "ERROR", "REFERENCE_ERROR"]
    }
    print(json.dumps(summary))
    for name, result in report["cases"].items():
        if result["status"] != "PASS":
            print(name, result.get("failures", result.get("reason", result.get("groups"))))


if __name__ == "__main__":
    main()
