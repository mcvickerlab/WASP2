"""Locked allele counting for one final bulk BAM per donor."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Literal

import polars as pl
import pysam

from .run_counting import run_count_variants

AnalysisUnit = Literal["snv", "feature"]
SCHEMA_VERSION = "wasp2.count-bundle.v1"
UTC = timezone.utc


@dataclass(frozen=True)
class DonorInput:
    donor_id: str
    vcf_sample: str
    bam: Path
    bam_index: Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _signature(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns


def _require_unchanged(path: Path, signature: tuple[int, int, int, int]) -> None:
    if _signature(path) != signature:
        raise RuntimeError(f"Input changed while cohort counting was running: {path}")


def _file_record(path: Path) -> tuple[dict[str, str | int], tuple[int, int, int, int]]:
    resolved = path.resolve(strict=True)
    before = _signature(resolved)
    digest = _sha256(resolved)
    _require_unchanged(resolved, before)
    stat = resolved.stat()
    return (
        {
            "path": str(resolved),
            "sha256": digest,
            "size_bytes": stat.st_size,
        },
        before,
    )


def _require_matching_hash(path: Path, expected_sha256: str) -> None:
    if _sha256(path) != expected_sha256:
        raise RuntimeError(f"Input content changed while cohort counting was running: {path}")


def _bam_index(bam: Path) -> Path:
    candidates = [
        Path(f"{bam}.bai"),
        bam.with_suffix(".bai"),
        Path(f"{bam}.csi"),
        bam.with_suffix(".csi"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve(strict=True)
    raise ValueError(f"Missing BAM index for {bam}; expected a BAI or CSI sidecar")


def _validate_bam_index(bam: Path, index: Path) -> None:
    try:
        with pysam.AlignmentFile(str(bam), "rb", index_filename=str(index)) as alignment:
            if not alignment.has_index() or not alignment.check_index():
                raise ValueError(f"BAM index is not usable for {bam}")
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot open indexed BAM: {bam}") from error


def _variant_index(variants: Path) -> Path:
    name = variants.name.lower()
    if not (name.endswith(".vcf.gz") or name.endswith(".vcf.bgz") or name.endswith(".bcf")):
        raise ValueError("Cohort counting requires an indexed VCF.GZ, VCF.BGZ, or BCF file")
    candidates = [Path(f"{variants}.csi"), Path(f"{variants}.tbi")]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve(strict=True)
    raise ValueError(f"Missing variant index for {variants}; expected .csi or .tbi sidecar")


def _vcf_samples(variants: Path) -> set[str]:
    try:
        with pysam.VariantFile(str(variants)) as variant_file:
            next(variant_file.fetch(), None)
            return set(variant_file.header.samples)
    except (OSError, ValueError) as error:
        raise ValueError(f"Cannot read cohort variant file: {variants}") from error


def _read_donor_manifest(path: Path, valid_vcf_samples: set[str]) -> list[DonorInput]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["donor_id", "vcf_sample", "bam"]:
            raise ValueError(
                "Donor manifest must have exactly three columns: donor_id, vcf_sample, bam"
            )
        rows = list(reader)
    if not rows:
        raise ValueError("Donor manifest contains no donors")

    donors: list[DonorInput] = []
    seen_donors: set[str] = set()
    seen_samples: set[str] = set()
    seen_bams: set[Path] = set()
    for line_number, row in enumerate(rows, start=2):
        donor_id = row["donor_id"].strip()
        vcf_sample = row["vcf_sample"].strip()
        bam_value = row["bam"].strip()
        if not donor_id or not vcf_sample or not bam_value:
            raise ValueError(f"Empty donor_id, vcf_sample, or bam on manifest line {line_number}")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", donor_id):
            raise ValueError(f"Invalid donor_id: {donor_id!r}")
        if any(character in vcf_sample for character in "\t\r\n"):
            raise ValueError(f"Invalid vcf_sample: {vcf_sample!r}")
        if donor_id in seen_donors:
            raise ValueError(f"Duplicate donor_id in donor manifest: {donor_id}")
        if vcf_sample in seen_samples:
            raise ValueError(f"Duplicate vcf_sample in donor manifest: {vcf_sample}")
        if vcf_sample not in valid_vcf_samples:
            raise ValueError(f"VCF sample not found in cohort variant file: {vcf_sample}")

        bam = Path(bam_value).expanduser()
        if not bam.is_absolute():
            bam = path.parent / bam
        bam = bam.resolve(strict=True)
        if bam.suffix.lower() != ".bam":
            raise ValueError(f"Bulk cohort inputs must be BAM files: {bam}")
        if bam in seen_bams:
            raise ValueError(f"Each donor must have a distinct BAM: {bam}")

        bam_index = _bam_index(bam)
        _validate_bam_index(bam, bam_index)
        donors.append(DonorInput(donor_id, vcf_sample, bam, bam_index))
        seen_donors.add(donor_id)
        seen_samples.add(vcf_sample)
        seen_bams.add(bam)
    return sorted(donors, key=lambda donor: donor.donor_id)


def _package_version() -> str:
    try:
        return version("wasp2")
    except PackageNotFoundError:
        return "unknown"


def _git_state() -> dict[str, str | bool | None]:
    repository = Path(__file__).resolve().parents[2]
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repository,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=repository,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return {"git_commit": None, "git_dirty": None}
    commit = revision.stdout.strip()
    return {
        "git_commit": commit if re.fullmatch(r"[0-9a-f]{40}", commit) else None,
        "git_dirty": bool(status.stdout.strip()),
    }


def _normalize_snv_rows(frame: pl.DataFrame, donor_id: str) -> pl.DataFrame:
    required = {"chrom", "pos", "ref", "alt", "ref_count", "alt_count", "other_count"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise RuntimeError(f"Allele counter omitted required columns for {donor_id}: {missing}")

    identity = ["chrom", "pos", "ref", "alt"]
    compared = ["ref_count", "alt_count", "other_count"]
    if "GT" in frame.columns:
        compared.append("GT")
    conflicts = frame.group_by(identity).agg(
        [pl.col(column).n_unique().alias(f"n_{column}") for column in compared]
    )
    conflict_filter = pl.any_horizontal([pl.col(f"n_{column}") > 1 for column in compared])
    if conflicts.filter(conflict_filter).height:
        row = conflicts.filter(conflict_filter).row(0, named=True)
        raise RuntimeError(
            f"Conflicting count rows for {donor_id} "
            f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"
        )

    columns = [
        column
        for column in [
            "chrom",
            "pos0",
            "pos",
            "ref",
            "alt",
            "GT",
            "ref_count",
            "alt_count",
            "other_count",
        ]
        if column in frame.columns
    ]
    return (
        frame.select(columns)
        .unique(subset=identity, keep="first", maintain_order=True)
        .with_columns(pl.lit(donor_id).alias("donor_id"))
        .select(["donor_id", *columns])
        .sort(["chrom", "pos", "ref", "alt"])
    )


def _normalize_feature_rows(frame: pl.DataFrame, donor_id: str) -> pl.DataFrame:
    required = {
        "chrom",
        "pos",
        "ref",
        "alt",
        "GT",
        "region",
        "ref_count",
        "alt_count",
        "other_count",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise RuntimeError(f"Feature counter omitted required columns for {donor_id}: {missing}")
    if frame.get_column("region").is_null().any() or (frame.get_column("region") == "").any():
        raise RuntimeError(f"Feature counter produced an empty feature identifier for {donor_id}")
    return (
        frame.with_columns(pl.lit(donor_id).alias("donor_id"))
        .select(["donor_id", *frame.columns])
        .sort(["chrom", "pos", "ref", "alt", "region"])
    )


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key in count bundle manifest: {key}")
        result[key] = value
    return result


def _validate_hash(value: Any, label: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"Invalid SHA-256 for {label}")
    return value


def _bundle_member(bundle_dir: Path, record: dict[str, Any], expected_name: str) -> Path:
    if record.get("path") != expected_name:
        raise ValueError(f"Count bundle output path must be {expected_name}")
    candidate = bundle_dir / expected_name
    if candidate.is_symlink() or not candidate.is_file():
        raise ValueError(f"Count bundle output must be a regular non-symlink file: {expected_name}")
    return candidate


def verify_count_bundle(
    manifest_file: str | Path,
    *,
    expected_manifest_sha256: str | None = None,
) -> dict[str, Any]:
    """Verify locked count outputs and return the parsed manifest."""
    manifest_input = Path(manifest_file).expanduser()
    if manifest_input.is_symlink() or not manifest_input.is_file():
        raise ValueError(
            f"Count bundle manifest must be a regular non-symlink file: {manifest_input}"
        )
    manifest_path = manifest_input.resolve(strict=True)
    if expected_manifest_sha256 is not None:
        expected = _validate_hash(expected_manifest_sha256, "count bundle manifest")
        if _sha256(manifest_path) != expected:
            raise ValueError("Count bundle manifest failed expected SHA-256 verification")
    try:
        manifest = json.loads(
            manifest_path.read_text(),
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except (OSError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"Cannot read count bundle manifest: {manifest_path}") from error
    required_top_level = {
        "schema_version",
        "created_at_utc",
        "software",
        "analysis_unit",
        "statistical_grouping",
        "configuration",
        "inputs",
        "outputs",
    }
    if set(manifest) != required_top_level:
        raise ValueError("Count bundle manifest has missing or unknown top-level fields")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Unsupported count bundle schema: {manifest.get('schema_version')!r}")
    if manifest.get("analysis_unit") not in {"snv", "feature"}:
        raise ValueError("Count bundle has an invalid analysis_unit")

    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != {"counts", "donors"}:
        raise ValueError("Count bundle must declare exactly counts and donors outputs")
    resolved: dict[str, Path] = {}
    for name, record in outputs.items():
        if not isinstance(record, dict):
            raise ValueError(f"Invalid output record for {name}")
        expected_name = "counts.tsv.gz" if name == "counts" else "donors.tsv"
        path = _bundle_member(manifest_path.parent, record, expected_name)
        expected_hash = _validate_hash(record.get("sha256"), name)
        if path.stat().st_size != record.get("size_bytes") or _sha256(path) != expected_hash:
            raise ValueError(f"Count bundle output failed verification: {name}")
        resolved[name] = path

    try:
        counts = pl.read_csv(resolved["counts"], separator="\t")
        donors = pl.read_csv(resolved["donors"], separator="\t")
    except (OSError, pl.exceptions.PolarsError) as error:
        raise ValueError("Count bundle tables cannot be fully parsed") from error
    if not counts.columns or counts.columns[0] != "donor_id" or counts.is_empty():
        raise ValueError("Locked counts must begin with a donor_id column")
    if donors.columns != ["donor_id", "vcf_sample", "bam"] or donors.is_empty():
        raise ValueError("Locked donors.tsv has an invalid schema")
    if sum(donors.null_count().row(0)):
        raise ValueError("Locked donors.tsv contains missing values")
    for column in ["donor_id", "vcf_sample", "bam"]:
        if donors.get_column(column).n_unique() != donors.height:
            raise ValueError(f"Locked donors.tsv contains duplicate {column} values")
    count_donors = set(counts.get_column("donor_id").unique().to_list())
    donor_ids = set(donors.get_column("donor_id").to_list())
    if count_donors != donor_ids:
        raise ValueError("Locked count donor IDs do not match donors.tsv")
    if outputs["counts"].get("rows") != counts.height:
        raise ValueError("Locked count row total does not match the manifest")
    if outputs["counts"].get("donors") != donors.height:
        raise ValueError("Locked count donor total does not match donors.tsv")
    if outputs["donors"].get("rows") != donors.height:
        raise ValueError("Locked donor row total does not match the manifest")

    donor_inputs = manifest.get("inputs", {}).get("donors")
    if not isinstance(donor_inputs, list):
        raise ValueError("Count bundle donor provenance is missing")
    provenance_ids = [record.get("donor_id") for record in donor_inputs]
    if provenance_ids != sorted(donor_ids):
        raise ValueError("Count bundle donor provenance does not match donors.tsv")
    return manifest


def _fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish_bundle(staging: Path, destination: Path) -> None:
    """Claim the destination and publish the manifest as the commit marker."""
    destination.mkdir()
    try:
        for name in ["counts.tsv.gz", "donors.tsv", "count_manifest.json"]:
            os.replace(staging / name, destination / name)
        _fsync_directory(destination)
        staging.rmdir()
    except BaseException:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def run_count_cohort(
    donor_manifest: str,
    variant_file: str,
    output_dir: str,
    *,
    unit: AnalysisUnit | Literal["peak"],
    region_file: str | None = None,
    use_region_names: bool = True,
) -> dict[str, str | int]:
    """Count one final BAM per donor and publish a locked cohort bundle."""
    normalized_unit: AnalysisUnit = "feature" if unit == "peak" else unit
    if normalized_unit not in {"snv", "feature"}:
        raise ValueError(f"Unknown analysis unit: {unit}")

    manifest_path = Path(donor_manifest).expanduser().resolve(strict=True)
    variants = Path(variant_file).expanduser().resolve(strict=True)
    variants_index = _variant_index(variants)
    regions = Path(region_file).expanduser().resolve(strict=True) if region_file else None
    if normalized_unit == "feature" and regions is None:
        raise ValueError("Feature counting requires --regions")

    destination = Path(os.path.abspath(Path(output_dir).expanduser()))
    if os.path.lexists(destination):
        raise FileExistsError(f"Refusing to overwrite locked output directory: {destination}")
    if not destination.parent.is_dir():
        raise ValueError(f"Output parent directory does not exist: {destination.parent}")

    donors = _read_donor_manifest(manifest_path, _vcf_samples(variants))
    source_paths = [manifest_path, variants, variants_index]
    if regions is not None:
        source_paths.append(regions)
    source_paths.extend(path for donor in donors for path in [donor.bam, donor.bam_index])

    records: dict[Path, dict[str, str | int]] = {}
    signatures: dict[Path, tuple[int, int, int, int]] = {}
    for path in source_paths:
        record, signature = _file_record(path)
        records[path] = record
        signatures[path] = signature

    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent))
    work = staging / "_work"
    work.mkdir()
    try:
        normalized_donors = pl.DataFrame(
            {
                "donor_id": [donor.donor_id for donor in donors],
                "vcf_sample": [donor.vcf_sample for donor in donors],
                "bam": [str(donor.bam) for donor in donors],
            }
        )
        donors_path = staging / "donors.tsv"
        normalized_donors.write_csv(donors_path, separator="\t")

        donor_records: list[dict[str, Any]] = []
        counts_path = staging / "counts.tsv.gz"
        output_schema: pl.Schema | None = None
        total_rows = 0
        with counts_path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                for donor_index, donor in enumerate(donors):
                    donor_dir = work / donor.donor_id
                    donor_dir.mkdir()
                    donor_counts = donor_dir / "counts.tsv"
                    run_count_variants(
                        bam_file=str(donor.bam),
                        variant_file=str(variants),
                        region_file=str(regions) if regions else None,
                        samples=[donor.vcf_sample],
                        use_region_names=use_region_names,
                        out_file=str(donor_counts),
                        temp_loc=str(donor_dir),
                        include_indels=False,
                        biallelic_only=True,
                        use_rust=True,
                    )
                    _require_unchanged(donor.bam, signatures[donor.bam])
                    _require_unchanged(donor.bam_index, signatures[donor.bam_index])

                    frame = pl.read_csv(donor_counts, separator="\t")
                    if frame.is_empty():
                        raise RuntimeError(
                            f"Allele counting produced no rows for donor {donor.donor_id}"
                        )
                    if normalized_unit == "snv":
                        frame = _normalize_snv_rows(frame, donor.donor_id)
                    else:
                        frame = _normalize_feature_rows(frame, donor.donor_id)
                    if output_schema is None:
                        output_schema = frame.schema
                    elif frame.schema != output_schema:
                        raise RuntimeError(
                            f"Count schema for {donor.donor_id} differs from prior donors"
                        )
                    frame.write_csv(
                        compressed,
                        include_header=donor_index == 0,
                        separator="\t",
                    )
                    total_rows += frame.height
                    donor_records.append(
                        {
                            "donor_id": donor.donor_id,
                            "vcf_sample": donor.vcf_sample,
                            "bam": records[donor.bam],
                            "bam_index": records[donor.bam_index],
                            "count_rows": frame.height,
                        }
                    )

        for path, signature in signatures.items():
            _require_unchanged(path, signature)
            expected_hash = records[path]["sha256"]
            if not isinstance(expected_hash, str):
                raise RuntimeError(f"Missing locked input hash for {path}")
            _require_matching_hash(path, expected_hash)

        counts_record, _ = _file_record(counts_path)
        donors_record, _ = _file_record(donors_path)
        counts_record["path"] = "counts.tsv.gz"
        donors_record["path"] = "donors.tsv"
        input_records: dict[str, Any] = {
            "donor_manifest": records[manifest_path],
            "variants": records[variants],
            "variants_index": records[variants_index],
            "donors": donor_records,
        }
        if regions is not None:
            input_records["regions"] = records[regions]
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": datetime.now(UTC).isoformat(),
            "software": {
                "package": "wasp2",
                "version": _package_version(),
                **_git_state(),
            },
            "analysis_unit": normalized_unit,
            "statistical_grouping": (
                "donor-plus-variant" if normalized_unit == "snv" else "donor-plus-feature"
            ),
            "configuration": {
                "one_final_bam_per_donor": True,
                "biallelic_snv_only": True,
                "include_indels": False,
                "use_region_names": use_region_names,
            },
            "inputs": input_records,
            "outputs": {
                "counts": {
                    **counts_record,
                    "rows": total_rows,
                    "donors": len(donors),
                },
                "donors": {
                    **donors_record,
                    "rows": len(donors),
                },
            },
        }
        manifest_path_out = staging / "count_manifest.json"
        manifest_path_out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        for path in [counts_path, donors_path, manifest_path_out]:
            _fsync_file(path)
        _fsync_directory(staging)
        verify_count_bundle(manifest_path_out)

        shutil.rmtree(work)
        _publish_bundle(staging, destination)
        return {
            "output_dir": str(destination),
            "counts": str(destination / "counts.tsv.gz"),
            "donor_manifest": str(destination / "donors.tsv"),
            "manifest": str(destination / "count_manifest.json"),
            "donors": len(donors),
            "rows": total_rows,
        }
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
