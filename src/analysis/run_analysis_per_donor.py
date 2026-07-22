"""Donor-local allelic imbalance orchestration."""

from __future__ import annotations

import csv
import gzip
import hashlib
import inspect
import json
import math
import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Literal, cast

import pandas as pd

try:
    from wasp2_rust import analyze_imbalance_run as rust_analyze_imbalance_run
    from wasp2_rust import fit_imbalance_dispersion as rust_fit_imbalance_dispersion
except ImportError:
    rust_analyze_imbalance_run = None
    rust_fit_imbalance_dispersion = None

AnalysisUnit = Literal["snv", "feature"]
AnalysisUnitInput = Literal["snv", "feature", "peak"]
DispersionScope = Literal["global", "per-donor"]
Model = Literal["single", "linear"]
FitMetadataValue = str | int | float | bool | None

UTC = timezone.utc
SCHEMA_VERSION = "wasp2.donor-local-analysis.v2"
COMMON_COLUMNS = {
    "donor_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "ref_count",
    "alt_count",
}
SNV_COLUMNS = ["donor_id", "chrom", "pos", "ref", "alt"]
RUST_RESULT_COLUMNS = [
    "region",
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
]
STATISTIC_COLUMNS = [
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
]


@dataclass(frozen=True)
class CountInput:
    """A verified count table and its optional locked-bundle manifest."""

    counts: Path
    counts_record: dict[str, str | int]
    counts_signature: tuple[int, int, int, int]
    manifest: Path | None
    manifest_record: dict[str, str | int] | None
    manifest_signature: tuple[int, int, int, int] | None
    manifest_data: dict[str, Any] | None

    @property
    def bundled(self) -> bool:
        return self.manifest is not None


@dataclass(frozen=True)
class DispersionFit:
    """Validated nuisance-parameter fit returned by Rust."""

    payload: Mapping[str, Any]
    fit_id: str
    method: Model
    rho: float | None
    linear_d1: float | None
    linear_d2: float | None
    linear_depth_center: float | None
    linear_depth_scale: float | None
    n_observations: int
    metadata: Mapping[str, FitMetadataValue]
    extended_contract: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _signature(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns


def _file_record(path: Path) -> tuple[dict[str, str | int], tuple[int, int, int, int]]:
    resolved = path.resolve(strict=True)
    before = _signature(resolved)
    digest = _sha256(resolved)
    if _signature(resolved) != before:
        raise RuntimeError(f"Input changed while it was being hashed: {resolved}")
    stat = resolved.stat()
    return (
        {
            "path": str(resolved),
            "sha256": digest,
            "size_bytes": stat.st_size,
            "mtime_utc": datetime.fromtimestamp(stat.st_mtime, UTC).isoformat(),
        },
        before,
    )


def _package_version() -> str:
    try:
        return version("wasp2")
    except PackageNotFoundError:
        return "unknown"


def _normalize_unit(unit: AnalysisUnitInput) -> AnalysisUnit:
    normalized = "feature" if unit == "peak" else unit
    if normalized not in {"snv", "feature"}:
        raise ValueError("unit must be 'snv' or 'feature' ('peak' is an alias)")
    return cast(AnalysisUnit, normalized)


def _validate_sha256(value: str, label: str) -> str:
    normalized = value.lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise ValueError(f"{label} must be a 64-character hexadecimal SHA-256 digest")
    return normalized


def _locate_manifest(source: Path) -> Path | None:
    if source.is_dir():
        manifest = source / "count_manifest.json"
        if not manifest.is_file():
            raise ValueError(f"Count bundle directory has no count_manifest.json: {source}")
        return manifest
    if source.name == "count_manifest.json":
        return source
    adjacent = source.parent / "count_manifest.json"
    return adjacent if adjacent.is_file() else None


def _resolve_count_input(
    count_file: str | Path,
    *,
    unit: AnalysisUnit,
    expected_manifest_sha256: str | None,
) -> CountInput:
    source = Path(count_file).expanduser()
    if not source.exists():
        raise FileNotFoundError(f"Count input does not exist: {source}")
    source = source.resolve(strict=True)
    manifest_path = _locate_manifest(source)

    if manifest_path is None:
        if source.is_dir():
            raise ValueError(f"Count input is not a regular file: {source}")
        if expected_manifest_sha256 is not None:
            raise ValueError(
                "--expected-manifest-sha256 requires a locked count bundle; "
                "standalone count tables have no manifest trust anchor"
            )
        counts_record, counts_signature = _file_record(source)
        return CountInput(
            counts=source,
            counts_record=counts_record,
            counts_signature=counts_signature,
            manifest=None,
            manifest_record=None,
            manifest_signature=None,
            manifest_data=None,
        )

    from counting.run_counting_cohort import verify_count_bundle

    expected = (
        _validate_sha256(expected_manifest_sha256, "expected manifest SHA-256")
        if expected_manifest_sha256 is not None
        else None
    )
    manifest_data = verify_count_bundle(
        manifest_path,
        expected_manifest_sha256=expected,
    )
    bundle_unit = manifest_data.get("analysis_unit")
    if bundle_unit != unit:
        raise ValueError(
            f"Requested analysis unit {unit!r} does not match locked count bundle unit "
            f"{bundle_unit!r}"
        )

    counts_path = (manifest_path.parent / "counts.tsv.gz").resolve(strict=True)
    if source.is_file() and source != manifest_path.resolve(strict=True) and source != counts_path:
        raise ValueError(
            f"Input is inside a count bundle but is not its locked counts table: {source}"
        )
    counts_record, counts_signature = _file_record(counts_path)
    manifest_record, manifest_signature = _file_record(manifest_path)
    return CountInput(
        counts=counts_path,
        counts_record=counts_record,
        counts_signature=counts_signature,
        manifest=manifest_path.resolve(strict=True),
        manifest_record=manifest_record,
        manifest_signature=manifest_signature,
        manifest_data=manifest_data,
    )


def _read_header(path: Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header:
        raise ValueError("Count table is empty or has no header")
    if any(not column for column in header):
        raise ValueError("Count table contains an empty column name")
    duplicate_columns = sorted({column for column in header if header.count(column) > 1})
    if duplicate_columns:
        raise ValueError(f"Count table contains duplicate columns: {duplicate_columns}")
    return header


def _coerce_nonnegative_integer(frame: pd.DataFrame, column: str, maximum: int) -> None:
    values = frame[column].astype(str)
    valid = values.str.fullmatch(r"[0-9]+")
    if not valid.all():
        raise ValueError(f"Column {column!r} contains invalid integer values")
    integers = values.map(int)
    if integers.gt(maximum).any():
        raise ValueError(f"Column {column!r} contains values larger than {maximum}")
    frame[column] = integers.astype("uint64")


def _validate_text_column(frame: pd.DataFrame, column: str) -> None:
    values = frame[column].astype(str)
    invalid = values.eq("") | values.str.contains(r"[\t\r\n]", regex=True)
    if invalid.any():
        raise ValueError(f"Column {column!r} contains missing, empty, or multiline values")
    frame[column] = values


def _normalize_genotypes(frame: pd.DataFrame) -> dict[str, int]:
    if "GT" not in frame.columns:
        raise ValueError("Phased feature analysis requires a GT column")
    gt = frame["GT"].astype(str)
    ref = frame["ref"]
    alt = frame["alt"]
    numeric_forward = gt.eq("0|1")
    numeric_reverse = gt.eq("1|0")
    allele_forward = gt.eq(ref + "|" + alt)
    allele_reverse = gt.eq(alt + "|" + ref)
    valid = numeric_forward | numeric_reverse | allele_forward | allele_reverse
    if not valid.all():
        bad = frame.loc[~valid, ["donor_id", "chrom", "pos", "ref", "alt", "GT"]].iloc[0]
        raise ValueError(
            f"Phased GT must be exact 0|1, 1|0, REF|ALT, or ALT|REF; found {bad.to_dict()}"
        )
    forward = numeric_forward | allele_forward
    frame["GT"] = forward.map({True: "0|1", False: "1|0"})
    return {
        "forward": int(forward.sum()),
        "reverse": int((~forward).sum()),
        "converted_from_alleles": int((allele_forward | allele_reverse).sum()),
    }


def _first_conflicting_group(
    frame: pd.DataFrame,
    identity: Sequence[str],
    compared: Sequence[str],
) -> tuple[Any, ...] | None:
    conflicts = (
        frame.groupby(list(identity), sort=False, dropna=False)[list(compared)]
        .nunique(dropna=False)
        .gt(1)
        .any(axis=1)
    )
    if not conflicts.any():
        return None
    key = conflicts.index[conflicts][0]
    return key if isinstance(key, tuple) else (key,)


def _load_counts(
    count_input: CountInput,
    *,
    unit: AnalysisUnit,
    region_col: str | None,
    phased: bool,
) -> tuple[pd.DataFrame, str | None, dict[str, Any]]:
    header = _read_header(count_input.counts)
    has_donor_id = "donor_id" in header
    has_sample = "sample" in header
    has_snv_id = "snv_id" in header
    if has_donor_id and has_sample:
        raise ValueError("Count table cannot contain both donor_id and legacy sample columns")
    if count_input.bundled and not has_donor_id:
        raise ValueError("Locked count bundles require the canonical donor_id column")
    if not has_donor_id and not has_sample:
        raise ValueError("Count table requires a donor_id column")

    selected_region = region_col
    if unit == "snv" and region_col is not None:
        raise ValueError(
            "SNV analysis tests variants independently and does not accept --region-col"
        )
    if unit == "feature":
        if count_input.bundled and region_col not in {None, "region"}:
            raise ValueError("Locked feature bundles require the canonical region column")
        selected_region = "region" if region_col is None and count_input.bundled else region_col
        if not selected_region:
            raise ValueError(
                "Standalone feature analysis requires --region-col; locked feature bundles use region"
            )
        if selected_region in COMMON_COLUMNS or selected_region in {"sample", "snv_id"}:
            raise ValueError(
                f"Feature region column conflicts with a reserved column: {selected_region}"
            )

    required = set(COMMON_COLUMNS)
    if not has_donor_id:
        required.remove("donor_id")
        required.add("sample")
    if selected_region is not None:
        required.add(selected_region)
    missing = sorted(required.difference(header))
    if missing:
        raise ValueError(f"Count table is missing required columns: {missing}")

    frame = pd.read_csv(
        count_input.counts,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
        low_memory=False,
    )
    if frame.empty:
        raise ValueError("Count table contains no observations")
    if has_sample:
        frame = frame.rename(columns={"sample": "donor_id"})

    text_columns = ["donor_id", "chrom", "ref", "alt"]
    if has_snv_id:
        text_columns.append("snv_id")
    if selected_region is not None:
        text_columns.append(selected_region)
    for column in text_columns:
        _validate_text_column(frame, column)
    for column in ["pos", "ref_count", "alt_count"]:
        _coerce_nonnegative_integer(frame, column, 2**32 - 1)
    if frame["pos"].eq(0).any():
        raise ValueError("SNV positions must be one-based")
    invalid_alleles = (
        frame["ref"].eq(frame["alt"])
        | frame["ref"].str.len().ne(1)
        | frame["alt"].str.len().ne(1)
        | frame["alt"].str.contains(",", regex=False)
    )
    if invalid_alleles.any():
        raise ValueError("Donor-local analysis requires biallelic SNV observations")

    site_identity = ["donor_id", "chrom", "pos"]
    allele_conflict = _first_conflicting_group(frame, site_identity, ["ref", "alt"])
    if allele_conflict is not None:
        donor_id, chrom, pos = allele_conflict
        raise ValueError(f"Conflicting alleles for donor/SNV {donor_id} {chrom}:{pos}")

    if has_snv_id:
        snv_id_conflict = _first_conflicting_group(
            frame,
            ["donor_id", "snv_id"],
            ["chrom", "pos", "ref", "alt"],
        )
        if snv_id_conflict is not None:
            donor_id, snv_id = snv_id_conflict
            raise ValueError(f"snv_id {snv_id!r} maps to conflicting variants for donor {donor_id}")
        identity_conflict = _first_conflicting_group(frame, SNV_COLUMNS, ["snv_id"])
        if identity_conflict is not None:
            donor_id, chrom, pos, ref, alt = identity_conflict
            raise ValueError(
                f"Multiple snv_id values map to donor/SNV {donor_id} {chrom}:{pos}:{ref}>{alt}"
            )
    else:
        frame["snv_id"] = (
            frame["chrom"]
            + ":"
            + frame["pos"].astype(str)
            + ":"
            + frame["ref"]
            + ":"
            + frame["alt"]
        )

    gt_qc = _normalize_genotypes(frame) if phased else None
    repeated_values = ["ref_count", "alt_count"]
    if phased:
        repeated_values.append("GT")
    count_conflict = _first_conflicting_group(frame, SNV_COLUMNS, repeated_values)
    if count_conflict is not None:
        donor_id, chrom, pos, ref, alt = count_conflict
        raise ValueError(
            f"Conflicting repeated rows for donor/SNV {donor_id} {chrom}:{pos}:{ref}>{alt}"
        )

    raw_rows = len(frame)
    identity = [*SNV_COLUMNS, selected_region] if selected_region is not None else SNV_COLUMNS
    frame = frame.drop_duplicates(identity, keep="first").copy()
    frame["pos"] = frame["pos"].astype("uint32")
    frame["ref_count"] = frame["ref_count"].astype("uint32")
    frame["alt_count"] = frame["alt_count"].astype("uint32")
    frame["_total_count"] = frame["ref_count"].astype("uint64") + frame["alt_count"].astype(
        "uint64"
    )
    return (
        frame,
        selected_region,
        {
            "raw_rows": raw_rows,
            "validated_rows": len(frame),
            "exact_duplicate_rows_removed": raw_rows - len(frame),
            "donors": int(frame["donor_id"].nunique()),
            "unique_donor_snv_rows": int(frame.drop_duplicates(SNV_COLUMNS).shape[0]),
            "phase": gt_qc,
            "legacy_sample_normalized": has_sample,
            "snv_id_source": "input" if has_snv_id else "synthesized",
        },
    )


def _artifact_paths(result_path: Path) -> dict[str, Path]:
    if result_path.suffix != ".tsv":
        raise ValueError("Donor-local output must use a .tsv filename")
    stem = result_path.stem
    return {
        "results": result_path,
        "dispersion": result_path.with_name(f"{stem}.dispersion.tsv"),
        "qc": result_path.with_name(f"{stem}.qc.tsv"),
        "provenance": result_path.with_name(f"{stem}.provenance.json"),
    }


def _write_tsv(frame: pd.DataFrame, path: Path, columns: Sequence[str]) -> None:
    frame.loc[:, list(columns)].to_csv(path, sep="\t", header=True, index=False)


def _finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise RuntimeError(f"Rust returned an invalid {label}")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Rust returned an invalid {label}") from error
    if not math.isfinite(result):
        raise RuntimeError(f"Rust returned a non-finite {label}")
    return result


def _fit_metadata(
    raw_fit: Mapping[str, Any], *, add_parameters_returned_status: bool
) -> dict[str, FitMetadataValue]:
    metadata: dict[str, FitMetadataValue] = {}
    for name, value in raw_fit.items():
        if name != "fit_status" and name != "optimizer" and not name.startswith("optimizer_"):
            continue
        if value is None:
            continue
        if isinstance(value, str):
            if not value or any(character in value for character in "\t\r\n"):
                raise RuntimeError(f"Rust returned an invalid {name}")
            metadata[name] = value
        elif isinstance(value, (bool, int)):
            metadata[name] = value
        elif isinstance(value, float):
            if not math.isfinite(value):
                raise RuntimeError(f"Rust returned a non-finite {name}")
            metadata[name] = value
        else:
            raise RuntimeError(f"Rust returned a non-scalar {name}")
    if add_parameters_returned_status and "fit_status" not in metadata:
        # The API currently proves that parameters were returned, but does not
        # expose optimizer convergence diagnostics.
        metadata["fit_status"] = "parameters_returned"
    return metadata


def _validate_fit(
    raw_fit: Mapping[str, Any],
    *,
    model: Model,
    expected_rows: int,
    allow_legacy_linear_mock: bool = False,
) -> DispersionFit:
    fit_id = raw_fit.get("fit_id")
    if not isinstance(fit_id, str) or not fit_id:
        raise RuntimeError("Rust dispersion fit did not return a non-empty fit_id")
    returned_method = raw_fit.get("method", raw_fit.get("model"))
    if returned_method != model:
        raise RuntimeError(
            f"Rust dispersion fit method {returned_method!r} does not match requested {model!r}"
        )
    try:
        n_observations = int(raw_fit["n_observations"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError("Rust dispersion fit did not return n_observations") from error
    if n_observations != expected_rows:
        raise RuntimeError(
            f"Rust dispersion fit used {n_observations} observations; expected {expected_rows}"
        )

    if model == "single":
        rho = _finite_float(raw_fit.get("rho"), "rho")
        linear_d1 = None
        linear_d2 = None
        linear_depth_center = None
        linear_depth_scale = None
        linear_values = (
            raw_fit.get("linear_d1"),
            raw_fit.get("linear_d2"),
            raw_fit.get("linear_depth_center"),
            raw_fit.get("linear_depth_scale"),
        )
        if any(value is not None for value in linear_values):
            raise RuntimeError("Rust single-dispersion fit returned linear parameters")
        extended_contract = True
    else:
        rho = None
        linear_d1 = _finite_float(raw_fit.get("linear_d1"), "linear_d1")
        linear_d2 = _finite_float(raw_fit.get("linear_d2"), "linear_d2")
        if raw_fit.get("rho") is not None:
            raise RuntimeError("Rust linear-dispersion fit returned scalar rho")
        raw_center = raw_fit.get("linear_depth_center")
        raw_scale = raw_fit.get("linear_depth_scale")
        if (raw_center is None) != (raw_scale is None):
            raise RuntimeError(
                "Rust linear-dispersion fit returned incomplete depth standardization"
            )
        if raw_center is None:
            if not allow_legacy_linear_mock:
                raise RuntimeError(
                    "Rust linear-dispersion fit omitted linear_depth_center and linear_depth_scale"
                )
            linear_depth_center = None
            linear_depth_scale = None
            extended_contract = False
        else:
            linear_depth_center = _finite_float(raw_center, "linear_depth_center")
            linear_depth_scale = _finite_float(raw_scale, "linear_depth_scale")
            if linear_depth_scale <= 0:
                raise RuntimeError("Rust returned a non-positive linear_depth_scale")
            extended_contract = True
    metadata = _fit_metadata(
        raw_fit,
        add_parameters_returned_status=extended_contract,
    )
    return DispersionFit(
        payload=dict(raw_fit),
        fit_id=fit_id,
        method=model,
        rho=rho,
        linear_d1=linear_d1,
        linear_d2=linear_d2,
        linear_depth_center=linear_depth_center,
        linear_depth_scale=linear_depth_scale,
        n_observations=n_observations,
        metadata=metadata,
        extended_contract=extended_contract,
    )


def _generated_fit_id(
    path: Path,
    fit: Mapping[str, Any],
    *,
    model: Model,
    min_count: int,
    pseudocount: int,
) -> str:
    def float_identity(value: Any) -> str | None:
        return None if value is None else _finite_float(value, "dispersion parameter").hex()

    identity = {
        "input_sha256": _sha256(path),
        "method": model,
        "min_count": min_count,
        "pseudocount": pseudocount,
        "rho": float_identity(fit.get("rho")),
        "linear_d1": float_identity(fit.get("linear_d1")),
        "linear_d2": float_identity(fit.get("linear_d2")),
        "linear_depth_center": float_identity(fit.get("linear_depth_center")),
        "linear_depth_scale": float_identity(fit.get("linear_depth_scale")),
        "n_observations": fit.get("n_observations"),
    }
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _fit_dispersion(
    frame: pd.DataFrame,
    path: Path,
    *,
    model: Model,
    min_count: int,
    pseudocount: int,
) -> DispersionFit:
    fitter = rust_fit_imbalance_dispersion
    if fitter is None:
        raise RuntimeError("Rust dispersion-fitting API is not available")
    fit_columns = [*SNV_COLUMNS, "ref_count", "alt_count"]
    _write_tsv(frame, path, fit_columns)
    base_kwargs: dict[str, Any] = {
        "min_count": min_count,
        "pseudocount": pseudocount,
        "method": model,
    }
    full_kwargs = {**base_kwargs, "phased": False, "region_col": "donor_id"}
    is_pyo3 = _is_known_pyo3_api(fitter)
    if is_pyo3:
        # donor_id is part of the fit identity. This is essential for a global
        # fit because the same genomic site in two donors is two observations.
        raw_fit = fitter(str(path), **full_kwargs)
    else:
        parameters = inspect.signature(fitter).parameters
        accepts_keywords = any(
            parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in parameters.values()
        )
        supports_full_call = accepts_keywords or {"phased", "region_col"}.issubset(parameters)
        raw_fit = fitter(str(path), **(full_kwargs if supports_full_call else base_kwargs))
    if not isinstance(raw_fit, Mapping):
        raise RuntimeError("Rust dispersion-fitting API returned an invalid payload")
    payload = dict(raw_fit)
    payload.setdefault(
        "fit_id",
        _generated_fit_id(
            path,
            payload,
            model=model,
            min_count=min_count,
            pseudocount=pseudocount,
        ),
    )
    return _validate_fit(
        payload,
        model=model,
        expected_rows=len(frame),
        allow_legacy_linear_mock=not is_pyo3,
    )


def _is_known_pyo3_api(callback: Any) -> bool:
    module = getattr(callback, "__module__", "")
    return module == "wasp2_rust" or module.startswith("wasp2_rust.")


def _parameter_values(source: Mapping[str, Any]) -> tuple[Any, Any, Any, Any, Any]:
    return (
        source.get("rho"),
        source.get("linear_d1"),
        source.get("linear_d2"),
        source.get("linear_depth_center"),
        source.get("linear_depth_scale"),
    )


def _parameter_identity(values: Sequence[Any]) -> tuple[str | None, ...]:
    return tuple(
        None if value is None else _finite_float(value, "dispersion parameter").hex()
        for value in values
    )


def _call_analyzer_with_fit(
    path: Path,
    *,
    fit: DispersionFit,
    model: Model,
    phased: bool,
    region_col: str | None,
    min_count: int,
    pseudocount: int,
) -> Mapping[str, Any]:
    analyzer = rust_analyze_imbalance_run
    if analyzer is None:
        raise RuntimeError("Rust donor-analysis API is not available")
    base_kwargs: dict[str, Any] = {
        "min_count": min_count,
        "pseudocount": pseudocount,
        "phased": phased,
        "region_col": region_col,
        "method": model,
    }
    explicit_kwargs = {
        **base_kwargs,
        "rho": fit.rho,
        "linear_d1": fit.linear_d1,
        "linear_d2": fit.linear_d2,
        "linear_depth_center": fit.linear_depth_center,
        "linear_depth_scale": fit.linear_depth_scale,
        "exact_snv": region_col is None,
    }
    if _is_known_pyo3_api(analyzer):
        raw_run = analyzer(str(path), **explicit_kwargs)
    else:
        parameters = inspect.signature(analyzer).parameters
        if "dispersion_fit" in parameters or not {
            "rho",
            "linear_d1",
            "linear_d2",
        }.intersection(parameters):
            run = analyzer(str(path), **base_kwargs, dispersion_fit=fit.payload)
            if not isinstance(run, Mapping):
                raise RuntimeError("Rust donor-analysis API returned an invalid payload")
            return run
        raw_run = analyzer(str(path), **explicit_kwargs)

    if not isinstance(raw_run, Mapping):
        raise RuntimeError("Rust donor-analysis API returned an invalid payload")
    if _parameter_identity(_parameter_values(raw_run)) != _parameter_identity(
        _parameter_values(fit.payload)
    ):
        raise RuntimeError("Rust donor-analysis API changed the fixed dispersion parameters")
    if raw_run.get("nuisance_source") not in {None, "fixed"}:
        raise RuntimeError("Rust donor-analysis API refitted the supplied dispersion parameters")
    run = {**raw_run, "dispersion_fit": fit.payload}
    if not isinstance(run, Mapping):
        raise RuntimeError("Rust donor-analysis API returned an invalid payload")
    return run


def _run_fit_values(run: Mapping[str, Any]) -> tuple[Any, Any, Any, Any, Any, Any]:
    nested = run.get("dispersion_fit")
    source = nested if isinstance(nested, Mapping) else run
    fit_id = source.get("fit_id", source.get("dispersion_fit_id"))
    return fit_id, *_parameter_values(source)


def _require_exact_fit_echo(run: Mapping[str, Any], fit: DispersionFit, donor_id: str) -> None:
    observed = _run_fit_values(run)
    expected = (
        fit.fit_id,
        fit.rho,
        fit.linear_d1,
        fit.linear_d2,
        fit.linear_depth_center,
        fit.linear_depth_scale,
    )
    if observed[0] != expected[0] or _parameter_identity(observed[1:]) != _parameter_identity(
        expected[1:]
    ):
        raise RuntimeError(
            f"Rust donor run for {donor_id} did not reuse the requested dispersion fit exactly: "
            f"expected {expected!r}, observed {observed!r}"
        )


def _validate_rust_results(run: Mapping[str, Any], donor_id: str) -> pd.DataFrame:
    raw_results = run.get("results")
    if not isinstance(raw_results, list) or not raw_results:
        raise RuntimeError(f"Rust returned no results for included donor {donor_id}")
    results = pd.DataFrame(raw_results)
    missing = sorted(set(RUST_RESULT_COLUMNS).difference(results.columns))
    if missing:
        raise RuntimeError(f"Rust results for {donor_id} are missing columns: {missing}")
    if results["region"].astype(str).duplicated().any():
        raise RuntimeError(f"Rust returned duplicate result identities for donor {donor_id}")
    return results.loc[:, RUST_RESULT_COLUMNS].copy()


def _prepare_analysis_table(
    donor_frame: pd.DataFrame,
    path: Path,
    *,
    unit: AnalysisUnit,
    region_col: str | None,
    phased: bool,
) -> tuple[str | None, pd.DataFrame]:
    columns = ["chrom", "pos", "ref", "alt"]
    rust_region_col = None
    if phased:
        columns.append("GT")
    if unit == "feature":
        if region_col is None:
            raise RuntimeError("Feature analysis has no normalized region column")
        columns.append(region_col)
        rust_region_col = region_col
    columns.extend(["ref_count", "alt_count"])
    _write_tsv(donor_frame, path, columns)
    return rust_region_col, donor_frame


def _format_donor_results(
    run: Mapping[str, Any],
    donor_frame: pd.DataFrame,
    *,
    donor_id: str,
    unit: AnalysisUnit,
    region_col: str | None,
) -> pd.DataFrame:
    results = _validate_rust_results(run, donor_id)
    try:
        n_observations = int(run["n_observations"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"Rust donor run for {donor_id} omitted n_observations") from error
    if n_observations != len(donor_frame):
        raise RuntimeError(
            f"Rust donor run for {donor_id} used {n_observations} rows; expected {len(donor_frame)}"
        )

    if unit == "snv":
        annotations = donor_frame.loc[
            :, ["snv_id", "chrom", "pos", "ref", "alt", "ref_count", "alt_count"]
        ].rename(
            columns={
                "ref_count": "expected_ref_count",
                "alt_count": "expected_alt_count",
            }
        )
        annotations["region"] = annotations["chrom"] + "_" + annotations["pos"].astype(str)
        results = annotations.merge(results, on="region", how="inner", validate="one_to_one")
        if len(results) != len(annotations):
            raise RuntimeError(f"Rust omitted eligible SNV results for donor {donor_id}")
        mismatch = (
            results["ref_count"].ne(results["expected_ref_count"])
            | results["alt_count"].ne(results["expected_alt_count"])
            | results["N"].ne(results["expected_ref_count"] + results["expected_alt_count"])
            | results["snp_count"].ne(1)
        )
        if mismatch.any():
            raise RuntimeError(f"Rust SNV counts differ from input for donor {donor_id}")
        leading_columns = ["donor_id", "snv_id", "chrom", "pos", "ref", "alt"]
        results = results.drop(columns=["region", "expected_ref_count", "expected_alt_count"])
    else:
        if region_col is None:
            raise RuntimeError("Feature analysis has no normalized region column")
        expected = (
            donor_frame.groupby(region_col, sort=False)
            .agg(
                expected_ref_count=("ref_count", "sum"),
                expected_alt_count=("alt_count", "sum"),
                expected_snp_count=("pos", "size"),
            )
            .reset_index()
            .rename(columns={region_col: "region"})
        )
        results = expected.merge(results, on="region", how="inner", validate="one_to_one")
        if len(results) != len(expected):
            raise RuntimeError(f"Rust omitted eligible feature results for donor {donor_id}")
        mismatch = (
            results["ref_count"].ne(results["expected_ref_count"])
            | results["alt_count"].ne(results["expected_alt_count"])
            | results["snp_count"].ne(results["expected_snp_count"])
        )
        if mismatch.any():
            raise RuntimeError(f"Rust feature aggregation differs from input for donor {donor_id}")
        results = results.drop(
            columns=["expected_ref_count", "expected_alt_count", "expected_snp_count"]
        ).rename(columns={"region": "feature_id"})
        leading_columns = ["donor_id", "feature_id"]

    results.insert(0, "donor_id", donor_id)
    results["significant_q05"] = results["fdr_pval"] < 0.05
    results["significant_q10"] = results["fdr_pval"] < 0.10
    return results.loc[
        :, [*leading_columns, *STATISTIC_COLUMNS, "significant_q05", "significant_q10"]
    ]


def _fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _same_inode(left: Path, right: Path) -> bool:
    try:
        left_stat = left.stat()
        right_stat = right.stat()
    except FileNotFoundError:
        return False
    return (left_stat.st_dev, left_stat.st_ino) == (right_stat.st_dev, right_stat.st_ino)


def _publish_no_replace(staged: Mapping[str, Path], destinations: Mapping[str, Path]) -> None:
    """Publish sidecars first and the result commit marker last without clobbering."""

    published: list[tuple[Path, Path]] = []
    try:
        for name in ["dispersion", "qc", "provenance", "results"]:
            source = staged[name]
            destination = destinations[name]
            os.link(source, destination, follow_symlinks=False)
            published.append((source, destination))
        _fsync_directory(destinations["results"].parent)
    except BaseException:
        for source, destination in reversed(published):
            if _same_inode(source, destination):
                destination.unlink()
        _fsync_directory(destinations["results"].parent)
        raise


def _require_unchanged(count_input: CountInput) -> None:
    if _signature(count_input.counts) != count_input.counts_signature:
        raise RuntimeError("Count input changed while donor-local analysis was running")
    if count_input.counts_record["sha256"] != _sha256(count_input.counts):
        raise RuntimeError("Count input content changed while donor-local analysis was running")
    if count_input.manifest is not None:
        if count_input.manifest_signature is None or count_input.manifest_record is None:
            raise RuntimeError("Internal error: bundle manifest lock is incomplete")
        if _signature(count_input.manifest) != count_input.manifest_signature:
            raise RuntimeError("Count manifest changed while donor-local analysis was running")
        if count_input.manifest_record["sha256"] != _sha256(count_input.manifest):
            raise RuntimeError(
                "Count manifest content changed while donor-local analysis was running"
            )


def run_per_donor_analysis(
    count_file: str | Path,
    out_file: str | Path,
    *,
    unit: AnalysisUnitInput,
    model: Model = "single",
    dispersion_scope: DispersionScope = "per-donor",
    region_col: str | None = None,
    phased: bool = False,
    min_count: int = 10,
    pseudocount: int = 0,
    min_donor_observations: int = 50,
    expected_manifest_sha256: str | None = None,
) -> dict[str, Path]:
    """Run independent donor tests with an explicit reusable dispersion fit."""

    normalized_unit = _normalize_unit(unit)
    if model not in {"single", "linear"}:
        raise ValueError("model must be 'single' or 'linear'")
    if dispersion_scope not in {"global", "per-donor"}:
        raise ValueError("dispersion_scope must be 'global' or 'per-donor'")
    if min_count < 0 or pseudocount < 0:
        raise ValueError("min_count and pseudocount must be nonnegative")
    if min_donor_observations < 1:
        raise ValueError("min_donor_observations must be positive")
    if normalized_unit == "snv" and phased:
        raise ValueError("Per-donor SNV analysis is unphased; --phased is not valid")

    count_input = _resolve_count_input(
        count_file,
        unit=normalized_unit,
        expected_manifest_sha256=expected_manifest_sha256,
    )
    frame, normalized_region_col, input_qc = _load_counts(
        count_input,
        unit=normalized_unit,
        region_col=region_col,
        phased=phased,
    )

    result_path = Path(os.path.abspath(Path(out_file).expanduser()))
    paths = _artifact_paths(result_path)
    if not result_path.parent.is_dir():
        raise ValueError(f"Output parent directory does not exist: {result_path.parent}")
    existing = [str(path) for path in paths.values() if os.path.lexists(path)]
    if existing:
        raise FileExistsError(f"Refusing to overwrite donor-local analysis artifacts: {existing}")
    if rust_fit_imbalance_dispersion is None or rust_analyze_imbalance_run is None:
        raise RuntimeError(
            "Rust donor-local analysis APIs are not available. Build the extension with "
            "`maturin develop --release` in the WASP2 environment."
        )

    eligible_rows = frame.loc[frame["_total_count"] >= min_count].copy()
    fit_rows = (
        eligible_rows.drop_duplicates(SNV_COLUMNS, keep="first")
        .sort_values(SNV_COLUMNS, kind="mergesort")
        .copy()
    )
    donor_qc = (
        frame.groupby("donor_id", sort=True)
        .agg(
            input_rows=("donor_id", "size"),
            unique_snv_rows=("pos", "size"),
        )
        .reset_index()
    )
    unique_snv_counts = (
        frame.drop_duplicates(SNV_COLUMNS)
        .groupby("donor_id", sort=True)
        .size()
        .rename("unique_snv_rows")
    )
    eligible_counts = (
        fit_rows.groupby("donor_id", sort=True).size().rename("eligible_snv_observations")
    )
    eligible_test_rows = (
        eligible_rows.groupby("donor_id", sort=True).size().rename("eligible_test_rows")
    )
    donor_qc = donor_qc.drop(columns="unique_snv_rows").join(unique_snv_counts, on="donor_id")
    donor_qc = donor_qc.join(eligible_counts, on="donor_id").join(eligible_test_rows, on="donor_id")
    donor_qc[["eligible_snv_observations", "eligible_test_rows"]] = (
        donor_qc[["eligible_snv_observations", "eligible_test_rows"]].fillna(0).astype(int)
    )
    donor_qc["analyzed"] = donor_qc["eligible_snv_observations"] >= min_donor_observations
    donor_qc["status"] = donor_qc["analyzed"].map(
        {True: "included", False: "excluded_min_observations"}
    )
    included_donors = donor_qc.loc[donor_qc["analyzed"], "donor_id"].tolist()
    if not included_donors:
        raise RuntimeError("No donors met the minimum eligible-SNV observation requirement")

    eligible_rows = eligible_rows.loc[eligible_rows["donor_id"].isin(included_donors)].copy()
    fit_rows = fit_rows.loc[fit_rows["donor_id"].isin(included_donors)].copy()
    staging = Path(tempfile.mkdtemp(prefix=f".{result_path.stem}.staging-", dir=result_path.parent))
    work = staging / "_work"
    work.mkdir()
    try:
        global_fit = None
        if dispersion_scope == "global":
            global_fit = _fit_dispersion(
                fit_rows,
                work / "global_fit.tsv",
                model=model,
                min_count=min_count,
                pseudocount=pseudocount,
            )

        result_frames: list[pd.DataFrame] = []
        dispersion_rows: list[dict[str, Any]] = []
        provenance_fits: dict[str, dict[str, Any]] = {}
        for donor_id in included_donors:
            donor_id = str(donor_id)
            donor_rows = eligible_rows.loc[eligible_rows["donor_id"] == donor_id].copy()
            donor_fit_rows = fit_rows.loc[fit_rows["donor_id"] == donor_id].copy()
            fit = global_fit or _fit_dispersion(
                donor_fit_rows,
                work / "donor_fit.tsv",
                model=model,
                min_count=min_count,
                pseudocount=pseudocount,
            )
            rust_region_col, donor_rows = _prepare_analysis_table(
                donor_rows,
                work / "donor_counts.tsv",
                unit=normalized_unit,
                region_col=normalized_region_col,
                phased=phased,
            )
            run = _call_analyzer_with_fit(
                work / "donor_counts.tsv",
                fit=fit,
                model=model,
                phased=phased,
                region_col=rust_region_col,
                min_count=min_count,
                pseudocount=pseudocount,
            )
            _require_exact_fit_echo(run, fit, donor_id)
            if bool(run.get("requested_phased")) != phased:
                raise RuntimeError(f"Rust donor run for {donor_id} changed requested_phased")
            if phased and not bool(run.get("effective_phased")):
                raise RuntimeError(f"Rust silently fell back to unphased analysis for {donor_id}")
            donor_results = _format_donor_results(
                run,
                donor_rows,
                donor_id=donor_id,
                unit=normalized_unit,
                region_col=normalized_region_col,
            )
            result_frames.append(donor_results)
            dispersion_row = {
                "donor_id": donor_id,
                "dispersion_scope": dispersion_scope,
                "fit_id": fit.fit_id,
                "model": model,
                "rho": fit.rho,
                "linear_d1": fit.linear_d1,
                "linear_d2": fit.linear_d2,
            }
            if fit.extended_contract:
                dispersion_row.update(
                    {
                        "linear_depth_center": fit.linear_depth_center,
                        "linear_depth_scale": fit.linear_depth_scale,
                    }
                )
            dispersion_row.update(fit.metadata)
            dispersion_row.update(
                {
                    "fit_observations": fit.n_observations,
                    "donor_eligible_snv_observations": len(donor_fit_rows),
                    "result_rows": len(donor_results),
                }
            )
            dispersion_rows.append(dispersion_row)

            if fit.extended_contract or fit.metadata:
                fit_record = {
                    "fit_id": fit.fit_id,
                    "model": model,
                    "rho": fit.rho,
                    "linear_d1": fit.linear_d1,
                    "linear_d2": fit.linear_d2,
                }
                if fit.extended_contract:
                    fit_record.update(
                        {
                            "linear_depth_center": fit.linear_depth_center,
                            "linear_depth_scale": fit.linear_depth_scale,
                        }
                    )
                fit_record.update(fit.metadata)
                fit_record["fit_observations"] = fit.n_observations
                existing_fit = provenance_fits.get(fit.fit_id)
                if existing_fit is None:
                    provenance_fits[fit.fit_id] = {**fit_record, "donor_ids": [donor_id]}
                else:
                    comparable = {
                        key: value for key, value in existing_fit.items() if key != "donor_ids"
                    }
                    if comparable != fit_record:
                        raise RuntimeError(f"Fit ID {fit.fit_id} mapped to inconsistent metadata")
                    existing_fit["donor_ids"].append(donor_id)

        _require_unchanged(count_input)
        sort_columns = (
            ["donor_id", "snv_id"]
            if normalized_unit == "snv"
            else [
                "donor_id",
                "feature_id",
            ]
        )
        results = pd.concat(result_frames, ignore_index=True).sort_values(
            sort_columns, kind="mergesort"
        )
        dispersion = pd.DataFrame(dispersion_rows).sort_values("donor_id", kind="mergesort")
        qc = donor_qc.sort_values("donor_id", kind="mergesort")

        staged_paths = {name: staging / path.name for name, path in paths.items()}
        results.to_csv(staged_paths["results"], sep="\t", header=True, index=False)
        dispersion.to_csv(staged_paths["dispersion"], sep="\t", header=True, index=False)
        qc.to_csv(staged_paths["qc"], sep="\t", header=True, index=False)
        for name in ["results", "dispersion", "qc"]:
            _fsync_file(staged_paths[name])

        output_records: dict[str, dict[str, str | int]] = {}
        for name in ["results", "dispersion", "qc"]:
            record, _ = _file_record(staged_paths[name])
            output_records[name] = {**record, "path": paths[name].name}
        provenance = {
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": datetime.now(UTC).isoformat(),
            "software": {"package": "wasp2", "version": _package_version()},
            "analysis": {
                "scope": "per-donor",
                "unit": normalized_unit,
                "model": model,
                "dispersion_scope": dispersion_scope,
                "region_col": normalized_region_col,
                "phased": phased,
                "multiple_testing": "Benjamini-Hochberg independently within each donor",
            },
            "parameters": {
                "min_count": min_count,
                "pseudocount": pseudocount,
                "min_donor_observations": min_donor_observations,
            },
            "inputs": {
                "counts": count_input.counts_record,
                "count_manifest": count_input.manifest_record,
                "count_bundle_schema": (
                    count_input.manifest_data.get("schema_version")
                    if count_input.manifest_data is not None
                    else None
                ),
                "qc": input_qc,
            },
            "donors": {
                "total": len(qc),
                "included": int(qc["analyzed"].sum()),
                "excluded": int((~qc["analyzed"]).sum()),
            },
            "outputs": output_records,
        }
        if provenance_fits:
            provenance["dispersion_fits"] = list(provenance_fits.values())
        staged_paths["provenance"].write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n"
        )
        _fsync_file(staged_paths["provenance"])
        _fsync_directory(staging)
        _require_unchanged(count_input)
        _publish_no_replace(staged_paths, paths)
        shutil.rmtree(staging, ignore_errors=True)
        return paths
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
