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
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.optimize import OptimizeResult
from scipy.sparse import csr_matrix
from scipy.special import expit

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
        # Preserve the archive's arithmetic as well as its likelihood. Reducing
        # pre-mixed probabilities with np.prod changes last-bit rounding.
        previous = 1
        for forward, reverse in zip(p1, p2, strict=True):
            forward_combined = previous * forward
            reverse_combined = previous * reverse
            previous = 0.5 * forward_combined + 0.5 * reverse_combined
        return float(first_ll - np.log(previous))

    module.opt_prob = scalar_prob
    module.opt_unphased_dp = aligned_dp


def fixed_bulk_reference(case: dict, ref: str = ARCHIVE_REF) -> pd.DataFrame:
    """Execute archived inference with supplied nuisance parameters, on this platform.

    Only nuisance fitting is replaced. The archived data preparation, effect
    optimizer, likelihood recurrence, LRT and BH functions remain the oracle.
    Linear predictions evaluate the supplied four-field API directly to avoid
    a coefficient-conversion roundtrip. Forwarding its ignored phase flag and
    scalar/array shape repairs are the same explicit reference adaptations.
    """
    modules = load_archive(ref)
    archive = modules["as_analysis"]
    fixed = case["fixed_dispersion"]
    pc = case["pseudocount"]
    if case["method"] == "single":
        original = archive.minimize_scalar

        def fixed_scalar_optimizer(fun, *args, **kwargs):
            if len(kwargs.get("args", ())) == 2:
                return OptimizeResult(x=fixed["rho"])
            return original(fun, *args, **kwargs)

        archive.minimize_scalar = fixed_scalar_optimizer
    else:
        adapt_linear_reference(modules)
        original_model = archive.linear_model
        original_parse = archive.parse_opt

        def fixed_linear_model(frame, region_col, phased=False):
            # frame already has the archive's pseudocounts, filtering and dedup.
            depth = frame["N"].to_numpy() - 2 * pc
            prediction = expit(
                fixed["linear_d1"]
                + ((depth - fixed["linear_depth_center"]) / fixed["linear_depth_scale"])
                * fixed["linear_d2"]
            )
            # The original model obtains coefficients from minimize, then calls
            # expit once to build this vector. Supply the fixed vector there.
            archive.minimize = lambda *args, **kwargs: OptimizeResult(
                x=[fixed["linear_d1"], fixed["linear_d2"]]
            )
            archive.expit = lambda _: prediction
            archive.parse_opt = lambda df, disp=None, phased=False: original_parse(
                df, disp=disp, phased=effective_phase
            )
            effective_phase = phased
            return original_model(frame, region_col, phased=phased)

        archive.linear_model = fixed_linear_model
    with contextlib.redirect_stdout(io.StringIO()):
        return archive.get_imbalance(
            pd.DataFrame(case["counts"]),
            method=case["method"],
            phased=case["phased"],
            pseudocount=pc,
            region_col="region",
        )


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


def prepare_bulk_reference_counts(
    counts: pd.DataFrame,
    *,
    min_count: int,
    pseudocount: int,
    effective_phased: bool,
    region_col: str | None = None,
) -> pd.DataFrame:
    """Reproduce archived filtering/deduplication while preserving the anchor order.

    Input counts are raw. Returned counts include pseudocounts and N, as used
    by the archived objective. Phase fallback must already have been resolved.
    """
    frame = counts.copy()
    if region_col is None:
        frame["region"] = frame["chrom"].astype(str) + "_" + frame["pos"].astype(str)
    elif region_col != "region":
        frame["region"] = frame[region_col]
    frame[["ref_count", "alt_count"]] = (
        frame[["ref_count", "alt_count"]].astype("int64") + pseudocount
    )
    frame["N"] = frame["ref_count"] + frame["alt_count"]
    frame = frame.loc[frame["N"] >= min_count + 2 * pseudocount].copy()
    columns = ["chrom", "pos", "ref_count", "alt_count", "N", "region"]
    if effective_phased:
        frame["GT"] = frame["GT"].str.split("|", regex=False).str[0].astype("uint8")
        columns.append("GT")
    return frame[columns].drop_duplicates()


def compare_bulk_equivalence(
    expected: pd.DataFrame,
    observed: pd.DataFrame,
    *,
    retained_counts: pd.DataFrame,
    effective_phased: bool,
    reference_log_likelihood: Callable[[str, float], float],
    dispersion_comparison: dict,
) -> dict:
    """Keep strict results; independently verify mirrored unphased bulk optima.

    The caller supplies the ordered, filtered/deduplicated objective observations,
    a Python reference objective, and the independently checked nuisance fit.
    This helper never edits results or applies to single-cell comparisons.
    """
    required = {"region", "ref_count", "alt_count", "N", "snp_count", *STAT_COLUMNS}
    if (
        expected.empty
        or not required.issubset(expected.columns)
        or not required.issubset(observed.columns)
    ):
        return {
            "strict_comparison": {"status": "FAIL", "reason": "missing required bulk outputs"},
            "statistical_equivalence": "FAIL",
            "failures": ["missing required bulk outputs"],
            "regions": [],
        }
    strict = compare_frames(expected, observed)
    failures = []
    evidence = []
    if dispersion_comparison.get("status") != "PASS":
        failures.append("dispersion comparison did not pass")
    if expected["region"].duplicated().any() or observed["region"].duplicated().any():
        failures.append("duplicate tested identities")
    if strict.get("reason") or any(item != "mu" for item in strict.get("failures", [])):
        failures.append("strict differences other than mu")
    if failures:
        return {
            "strict_comparison": strict,
            "statistical_equivalence": "FAIL",
            "failures": failures,
            "regions": evidence,
        }
    actual = observed.set_index("region")
    for row in expected.itertuples(index=False):
        identity = row.region
        python_mu = float(row.mu)
        rust_mu = float(actual.loc[identity, "mu"])
        finite = np.isfinite(python_mu) and np.isfinite(rust_mu)
        valid = finite and 0 <= python_mu <= 1 and 0 <= rust_mu <= 1
        same = valid and np.isclose(python_mu, rust_mu, rtol=5e-5, atol=5e-5)
        columns = {"region", "ref_count", "alt_count", "N"}
        rows = (
            retained_counts.loc[retained_counts["region"] == identity]
            if columns.issubset(retained_counts.columns)
            else pd.DataFrame()
        )
        symmetric = bool(
            not effective_phased
            and len(rows) > 1
            and len(rows) == row.snp_count
            and rows.iloc[0]["ref_count"] == rows.iloc[0]["alt_count"]
            and np.isfinite(rows[["ref_count", "alt_count", "N"]].to_numpy()).all()
            and (rows[["ref_count", "alt_count", "N"]].to_numpy() >= 0).all()
            and (rows["N"] == rows["ref_count"] + rows["alt_count"]).all()
        )
        item = {
            "region": identity,
            "python_mu": python_mu if finite else None,
            "rust_mu": rust_mu if finite else None,
            "symmetric_unphased_anchor": symmetric,
            "retained_observations": len(rows),
            "strict_mu_match": bool(same),
        }
        if len(rows) and np.isfinite(rows[["ref_count", "alt_count", "N"]].to_numpy()).all():
            item["first_retained_counts"] = {
                "ref_count": int(rows.iloc[0]["ref_count"]),
                "alt_count": int(rows.iloc[0]["alt_count"]),
                "N": int(rows.iloc[0]["N"]),
            }
        if same:
            item["status"] = "STRICT_MATCH"
        elif not valid or not symmetric:
            item["status"] = "FAIL"
            failures.append(f"{identity}: mu mismatch without a valid symmetric anchor")
        elif not np.isclose(rust_mu, 1 - python_mu, rtol=5e-5, atol=5e-5):
            item["status"] = "FAIL"
            failures.append(f"{identity}: proportions are not complementary")
        else:
            try:
                values = np.array(
                    [
                        reference_log_likelihood(identity, mu)
                        for mu in [python_mu, rust_mu, 1 - python_mu]
                    ],
                    dtype=float,
                )
                reported = np.array([row.alt_ll, actual.loc[identity, "alt_ll"]], dtype=float)
                equal = (
                    np.isfinite(values).all()
                    and np.allclose(values, values[0], rtol=5e-5, atol=5e-5)
                    and np.allclose(values[:2], reported, rtol=5e-5, atol=5e-5)
                )
                item["reference_log_likelihoods"] = (
                    values.tolist() if np.isfinite(values).all() else None
                )
                item["reference_evaluation_order"] = ["python_mu", "rust_mu", "1-python_mu"]
                item["status"] = "SYMMETRIC_EQUIVALENT" if equal else "FAIL"
                if not equal:
                    failures.append(f"{identity}: reference likelihoods do not agree")
            except Exception as exc:
                item["status"] = "FAIL"
                failures.append(f"{identity}: reference evaluation failed: {exc}")
        evidence.append(item)
    return {
        "strict_comparison": strict,
        "statistical_equivalence": "FAIL" if failures else "PASS",
        "failures": failures,
        "regions": evidence,
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
