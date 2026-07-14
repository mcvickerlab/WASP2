import hashlib
import json
import os
from pathlib import Path

import polars as pl
import pysam
import pytest
from typer.testing import CliRunner

from counting import count_alleles
from counting import run_counting_cohort as cohort
from counting.__main__ import app

SHARED_DATA = Path(__file__).parents[1] / "shared_data"


def _indexed_vcf(tmp_path: Path) -> Path:
    source = tmp_path / "variants.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tvcf_a\tvcf_b\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0\n"
    )
    compressed = tmp_path / "variants.vcf.gz"
    pysam.tabix_compress(str(source), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    variants = _indexed_vcf(tmp_path)
    regions = tmp_path / "peaks.bed"
    regions.write_text("chr1\t0\t100\tpeak1\n")
    rows = []
    for donor_id, vcf_sample in [("donor_b", "vcf_b"), ("donor_a", "vcf_a")]:
        bam = tmp_path / f"{donor_id}.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header) as alignment:
            read = pysam.AlignedSegment()
            read.query_name = donor_id
            read.query_sequence = "AAAAAAAAAA"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 9
            read.mapping_quality = 60
            read.cigar = ((0, 10),)
            read.query_qualities = pysam.qualitystring_to_array("FFFFFFFFFF")
            alignment.write(read)
        pysam.index(str(bam))
        rows.append(f"{donor_id}\t{vcf_sample}\t{bam.name}\n")
    manifest = tmp_path / "donors.input.tsv"
    manifest.write_text("donor_id\tvcf_sample\tbam\n" + "".join(rows))
    return manifest, variants, regions, tmp_path / "locked_counts"


def _fake_count(**kwargs) -> None:
    frame = pl.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos0": [9, 9],
            "pos": [10, 10],
            "ref": ["A", "A"],
            "alt": ["G", "G"],
            "GT": ["0|1", "0|1"],
            "region": ["peak1", "peak_overlap"],
            "ref_count": [8, 8],
            "alt_count": [4, 4],
            "other_count": [0, 0],
        }
    )
    frame.write_csv(kwargs["out_file"], separator="\t")


@pytest.mark.unit
def test_snv_bundle_locks_explicit_donor_mapping(tmp_path, monkeypatch):
    manifest, variants, regions, output = _inputs(tmp_path)
    calls = []

    def record_count(**kwargs) -> None:
        calls.append(kwargs)
        _fake_count(**kwargs)

    monkeypatch.setattr(cohort, "run_count_variants", record_count)
    result = cohort.run_count_cohort(
        str(manifest),
        str(variants),
        str(output),
        unit="snv",
        region_file=str(regions),
    )

    counts = pl.read_csv(result["counts"], separator="\t")
    assert counts.height == 2
    assert counts.columns[0] == "donor_id"
    assert counts.get_column("donor_id").to_list() == ["donor_a", "donor_b"]
    assert "region" not in counts.columns
    assert [call["samples"] for call in calls] == [["vcf_a"], ["vcf_b"]]

    donors = pl.read_csv(result["donor_manifest"], separator="\t")
    assert donors.columns == ["donor_id", "vcf_sample", "bam"]
    assert donors.get_column("donor_id").to_list() == ["donor_a", "donor_b"]

    locked = cohort.verify_count_bundle(result["manifest"])
    assert locked["schema_version"] == cohort.SCHEMA_VERSION
    assert locked["analysis_unit"] == "snv"
    assert locked["statistical_grouping"] == "donor-plus-variant"
    assert locked["configuration"]["biallelic_snv_only"] is True
    assert locked["configuration"]["include_indels"] is False
    assert locked["outputs"]["counts"]["sha256"]
    assert locked["inputs"]["variants_index"]["sha256"]
    assert locked["inputs"]["donors"][0]["bam_index"]["sha256"]
    assert sorted(path.name for path in output.iterdir()) == [
        "count_manifest.json",
        "counts.tsv.gz",
        "donors.tsv",
    ]


@pytest.mark.unit
@pytest.mark.parametrize("unit", ["feature", "peak"])
def test_feature_bundle_preserves_feature_rows_and_peak_alias(tmp_path, monkeypatch, unit):
    manifest, variants, regions, output = _inputs(tmp_path)
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    result = cohort.run_count_cohort(
        str(manifest),
        str(variants),
        str(output),
        unit=unit,
        region_file=str(regions),
    )

    counts = pl.read_csv(result["counts"], separator="\t")
    assert counts.height == 4
    assert counts.get_column("region").unique().sort().to_list() == ["peak1", "peak_overlap"]
    locked = cohort.verify_count_bundle(result["manifest"])
    assert locked["analysis_unit"] == "feature"
    assert locked["statistical_grouping"] == "donor-plus-feature"


@pytest.mark.unit
def test_feature_bundle_requires_regions(tmp_path):
    manifest, variants, _, output = _inputs(tmp_path)
    with pytest.raises(ValueError, match="requires --regions"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="feature")


@pytest.mark.unit
def test_manifest_requires_explicit_unique_mapping(tmp_path):
    manifest, variants, _, output = _inputs(tmp_path)
    bam = tmp_path / "donor_a.bam"
    manifest.write_text(
        f"donor_id\tvcf_sample\tbam\ndonor_a\tvcf_a\t{bam}\ndonor_a\tvcf_b\t{bam}\n"
    )
    with pytest.raises(ValueError, match="Duplicate donor_id"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")

    manifest.write_text("donor_id\tvcf_sample\tbam\ndonor_a\tnot_in_vcf\tdonor_a.bam\n")
    with pytest.raises(ValueError, match="VCF sample not found"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")


@pytest.mark.unit
def test_bundle_verifier_rejects_tampering(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    result = cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    Path(result["counts"]).write_bytes(Path(result["counts"]).read_bytes() + b"tampered")
    with pytest.raises(ValueError, match="failed verification: counts"):
        cohort.verify_count_bundle(result["manifest"])


@pytest.mark.unit
def test_failure_removes_staging_and_detects_input_mutation(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)

    def mutate_input(**kwargs) -> None:
        Path(kwargs["bam_file"]).write_bytes(b"changed")
        _fake_count(**kwargs)

    monkeypatch.setattr(cohort, "run_count_variants", mutate_input)
    with pytest.raises(RuntimeError, match="Input changed"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    assert not output.exists()
    assert not list(tmp_path.glob(".locked_counts.staging-*"))


@pytest.mark.unit
def test_cli_routes_locked_contract(monkeypatch):
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return {"donors": 2, "rows": 4, "output_dir": "/tmp/out"}

    monkeypatch.setattr("counting.__main__.run_count_cohort", fake_run)
    result = CliRunner().invoke(
        app,
        [
            "count-cohort",
            "donors.tsv",
            "variants.vcf.gz",
            "counts",
            "--unit",
            "peak",
            "--regions",
            "peaks.bed",
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["unit"] == "peak"
    assert captured["region_file"] == "peaks.bed"


@pytest.mark.unit
def test_bundle_refuses_overwrite(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    output.mkdir()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")


@pytest.mark.unit
def test_manifest_json_uses_relative_bundle_members(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    result = cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    locked = json.loads(Path(result["manifest"]).read_text())
    assert locked["outputs"]["counts"]["path"] == "counts.tsv.gz"
    assert locked["outputs"]["donors"]["path"] == "donors.tsv"


@pytest.mark.unit
def test_expected_manifest_digest_is_an_external_trust_anchor(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    result = cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    digest = hashlib.sha256(Path(result["manifest"]).read_bytes()).hexdigest()
    cohort.verify_count_bundle(result["manifest"], expected_manifest_sha256=digest)
    with pytest.raises(ValueError, match="failed expected SHA-256"):
        cohort.verify_count_bundle(result["manifest"], expected_manifest_sha256="0" * 64)


@pytest.mark.unit
def test_same_size_input_mutation_with_restored_mtime_is_detected(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    target = tmp_path / "donor_a.bam"

    def mutate_without_stat_change(**kwargs) -> None:
        if Path(kwargs["bam_file"]) == target:
            stat = target.stat()
            content = bytearray(target.read_bytes())
            content[-1] ^= 1
            target.write_bytes(content)
            os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        _fake_count(**kwargs)

    monkeypatch.setattr(cohort, "run_count_variants", mutate_without_stat_change)
    with pytest.raises(RuntimeError, match="Input content changed"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    assert not output.exists()


@pytest.mark.unit
def test_atomic_destination_claim_never_replaces_racing_directory(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    calls = 0

    def create_destination_during_count(**kwargs) -> None:
        nonlocal calls
        calls += 1
        _fake_count(**kwargs)
        if calls == 2:
            output.mkdir()

    monkeypatch.setattr(cohort, "run_count_variants", create_destination_during_count)
    with pytest.raises(FileExistsError):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    assert output.is_dir()
    assert not list(tmp_path.glob(".locked_counts.staging-*"))


@pytest.mark.unit
def test_broken_symlink_destination_is_treated_as_existing(tmp_path):
    manifest, variants, _, output = _inputs(tmp_path)
    output.symlink_to(tmp_path / "missing-target")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")


@pytest.mark.unit
def test_count_outputs_are_byte_deterministic(tmp_path, monkeypatch):
    manifest, variants, _, first = _inputs(tmp_path)
    second = tmp_path / "locked_counts_second"
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    first_result = cohort.run_count_cohort(str(manifest), str(variants), str(first), unit="snv")
    second_result = cohort.run_count_cohort(str(manifest), str(variants), str(second), unit="snv")
    assert Path(first_result["counts"]).read_bytes() == Path(second_result["counts"]).read_bytes()
    assert (
        Path(first_result["donor_manifest"]).read_bytes()
        == Path(second_result["donor_manifest"]).read_bytes()
    )


@pytest.mark.unit
def test_verifier_rejects_noncanonical_output_paths(tmp_path, monkeypatch):
    manifest, variants, _, output = _inputs(tmp_path)
    monkeypatch.setattr(cohort, "run_count_variants", _fake_count)
    result = cohort.run_count_cohort(str(manifest), str(variants), str(output), unit="snv")
    manifest_path = Path(result["manifest"])
    locked = json.loads(manifest_path.read_text())
    locked["outputs"]["counts"]["path"] = "./counts.tsv.gz"
    manifest_path.write_text(json.dumps(locked, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="must be counts.tsv.gz"):
        cohort.verify_count_bundle(manifest_path)


@pytest.mark.integration
@pytest.mark.parametrize("unit", ["snv", "feature"])
def test_real_counting_pipeline_publishes_verified_bundle(tmp_path, unit):
    if not count_alleles.RUST_AVAILABLE:
        pytest.skip("Rust counting extension is not installed")
    donor_manifest = tmp_path / "donors.tsv"
    donor_manifest.write_text(
        "donor_id\tvcf_sample\tbam\n"
        f"donor_b\tsample2\t{SHARED_DATA / 'sample2.bam'}\n"
        f"donor_a\tsample1\t{SHARED_DATA / 'sample1.bam'}\n"
    )
    output = tmp_path / f"real_{unit}_counts"
    result = cohort.run_count_cohort(
        str(donor_manifest),
        str(SHARED_DATA / "variants.vcf.gz"),
        str(output),
        unit=unit,
        region_file=str(SHARED_DATA / "regions.bed"),
    )
    locked = cohort.verify_count_bundle(result["manifest"])
    counts = pl.read_csv(result["counts"], separator="\t")
    assert counts.get_column("donor_id").unique().sort().to_list() == ["donor_a", "donor_b"]
    assert locked["outputs"]["counts"]["rows"] == counts.height
    assert ("region" in counts.columns) is (unit == "feature")
