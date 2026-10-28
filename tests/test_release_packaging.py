"""Publication guards must fail before changing files or uploading artifacts."""

import hashlib
import importlib.util
import io
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from urllib.error import HTTPError

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_export_verification_process_exits(tmp_path):
    """Export assertions must also finish normal threaded-library teardown."""
    pytest.importorskip("pyarrow")
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    sparse = pytest.importorskip("scipy.sparse")
    data = ad.AnnData(
        sparse.csr_matrix([[1], [1]]),
        obs=pd.DataFrame(index=["cellA-1", "cellB-1"]),
    )
    data.layers["ref"] = sparse.csr_matrix(np.array([[1], [0]]))
    data.layers["alt"] = sparse.csr_matrix(np.array([[0], [1]]))
    data.uns["feature"] = {"region": ["peak1"]}
    data.write_h5ad(tmp_path / "cells.h5ad")
    pd.DataFrame({"region": ["peak1"], "ref_count": [1], "pval": [0.5]}).to_csv(
        tmp_path / "phased.tsv", sep="\t", index=False
    )
    script = (ROOT / "scripts/container_smoke_test.sh").read_text()
    verifier = script.split("<<'PYSC'\n", 1)[1].split("\nPYSC", 1)[0]
    for _ in range(20):
        completed = subprocess.run(
            [sys.executable, "-u", "-", str(tmp_path)],
            input=verifier,
            text=True,
            capture_output=True,
            timeout=30,
            env={**os.environ, "WASP2_TEST_EXPORTS": "1"},
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert "PASS: Zarr and Parquet export round trips" in completed.stdout
        assert "PASS: exact single-cell counts" in completed.stdout


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("version", ["1.5.0banana", "1.5.0..", "1.5.0-rc.1", "01.5.0"])
def test_invalid_version_does_not_write(monkeypatch, version):
    release = load_script("release_version")
    monkeypatch.setattr(release, "write", lambda *args: pytest.fail("invalid version wrote files"))
    with pytest.raises(RuntimeError, match="invalid release version"):
        release.set_version(version, "2026-10-16")


def test_invalid_release_date_does_not_write(monkeypatch):
    release = load_script("release_version")
    monkeypatch.setattr(release, "write", lambda *args: pytest.fail("invalid date wrote files"))
    with pytest.raises(ValueError):
        release.set_version("1.5.0", "2026-02-30")


def test_publication_requires_matching_dated_metadata(monkeypatch):
    release = load_script("release_version")
    monkeypatch.setattr(release, "cargo_version", lambda: "1.5.0")
    files = {"CITATION.cff": 'version: "1.5.0"\n', "CHANGELOG.md": "## [Unreleased]\n"}
    monkeypatch.setattr(release, "read", files.__getitem__)
    assert release.check_release_metadata()
    files["CITATION.cff"] += 'date-released: "2026-10-16"\n'
    assert release.check_release_metadata()
    files["CHANGELOG.md"] = "## [1.5.0] - 2026-10-16\n"
    assert release.check_release_metadata() == []


@pytest.mark.parametrize("published", [False, True])
def test_pypi_resume_only_stages_missing_identical_artifacts(tmp_path, monkeypatch, published):
    upload = load_script("prepare_pypi_upload")
    packages = tmp_path / "packages"
    packages.mkdir()
    (packages / "wasp2-1.5.0.tar.gz").write_bytes(b"source")
    (packages / "wasp2-1.5.0.whl").write_bytes(b"wheel")
    if published:
        response = {
            "urls": [
                {
                    "filename": "wasp2-1.5.0.tar.gz",
                    "digests": {"sha256": hashlib.sha256(b"source").hexdigest()},
                }
            ]
        }
        monkeypatch.setattr(
            upload.urllib.request,
            "urlopen",
            lambda *a, **k: io.BytesIO(json.dumps(response).encode()),
        )
    else:

        def absent(*args, **kwargs):
            raise HTTPError("https://pypi.org", 404, "missing", {}, None)

        monkeypatch.setattr(upload.urllib.request, "urlopen", absent)
    destination = tmp_path / "upload"
    assert upload.stage(packages, destination, "1.5.0") == (1 if published else 2)
    assert (destination / "wasp2-1.5.0.whl").read_bytes() == b"wheel"
    assert (destination / "wasp2-1.5.0.tar.gz").exists() != published


def test_pypi_resume_refuses_different_bytes_before_staging(tmp_path, monkeypatch):
    upload = load_script("prepare_pypi_upload")
    packages = tmp_path / "packages"
    packages.mkdir()
    (packages / "wasp2.whl").write_bytes(b"different")
    response = {"urls": [{"filename": "wasp2.whl", "digests": {"sha256": "0" * 64}}]}
    monkeypatch.setattr(
        upload.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(json.dumps(response).encode())
    )
    destination = tmp_path / "upload"
    with pytest.raises(ValueError, match="Published artifact differs"):
        upload.stage(packages, destination, "1.5.0")
    assert not destination.exists()


def test_sbom_cannot_pass_with_an_empty_inventory(tmp_path):
    inventory = load_script("check_sbom")
    path = tmp_path / "sbom.json"
    path.write_text('{"components": []}')
    with pytest.raises(ValueError, match="missing SBOM components"):
        inventory.check(path, ["wasp2", "rust-htslib"])
    path.write_text('{"components": [{"name": "wasp2"}, {"name": "rust-htslib"}]}')
    inventory.check(path, ["wasp2", "rust-htslib"])


def native_manifest(tmp_path):
    manifest = tmp_path / "native-sources.json"
    manifest.write_text(
        json.dumps(
            {
                "components": [
                    {
                        "name": name,
                        "version": "1.25",
                        "sha256": "a" * 64,
                        "url": f"https://example.org/{name}-1.25.tar.bz2",
                    }
                    for name in ["htslib", "samtools", "bcftools"]
                ]
            }
        )
    )
    return manifest


def test_native_manifest_rejects_invalid_hash_before_importing(tmp_path, monkeypatch):
    native = load_script("native_sbom")
    manifest = native_manifest(tmp_path)
    contents = json.loads(manifest.read_text())
    contents["components"][0]["sha256"] = "not-a-source-hash"
    manifest.write_text(json.dumps(contents))
    monkeypatch.setattr(
        native, "_embedded_components", lambda: pytest.fail("invalid source imported")
    )
    with pytest.raises(ValueError, match="invalid source SHA-256"):
        native.collect(manifest)


def test_native_merge_preserves_scanner_inventory_and_independent_versions(tmp_path, monkeypatch):
    native = load_script("native_sbom")
    monkeypatch.setattr(
        native,
        "_embedded_components",
        lambda: [
            native._component("htslib", version, origin, "hts_version()")
            for origin, version in [("wasp2-rust", "1.19.1"), ("pysam", "1.24"), ("cyvcf2", "1.25")]
        ],
    )
    collected = native.collect(native_manifest(tmp_path))
    sources = [
        component for component in collected["components"] if "externalReferences" in component
    ]
    assert len(sources) == 3
    assert all(
        component["externalReferences"][0]["hashes"] == [{"alg": "SHA-256", "content": "a" * 64}]
        for component in sources
    )
    source_path = tmp_path / "native.json"
    source_path.write_text(json.dumps(collected))
    sbom = tmp_path / "sbom.json"
    original = {
        "bomFormat": "CycloneDX",
        "specVersion": "1.6",
        "metadata": {"tools": {"components": [{"name": "syft"}]}},
        "components": [{"type": "library", "name": "libc", "bom-ref": "scanner:libc"}],
        "dependencies": [{"ref": "scanner:libc", "dependsOn": []}],
    }
    sbom.write_text(json.dumps(original))
    native.merge(sbom, source_path)
    first = sbom.read_bytes()
    native.merge(sbom, source_path)
    assert sbom.read_bytes() == first
    result = json.loads(first)
    assert result["components"][0] == original["components"][0]
    assert result["metadata"] == original["metadata"]
    assert result["dependencies"] == original["dependencies"]
    copies = [component for component in result["components"] if component["name"] == "htslib"]
    assert len(copies) == 4
    assert len({component["bom-ref"] for component in copies}) == 4
    assert {component["version"] for component in copies} == {"1.19.1", "1.24", "1.25"}


@pytest.mark.parametrize("invalid", ["bad-json", "bad-hash", "conflicting-reference"])
def test_native_merge_rejects_invalid_input_without_overwriting(tmp_path, invalid):
    native = load_script("native_sbom")
    sbom = tmp_path / "sbom.json"
    original = native._component("htslib", "1.25", "pysam", "hts_version()")
    sbom.write_text(json.dumps({"bomFormat": "CycloneDX", "components": [original]}))
    before = sbom.read_bytes()
    additions = tmp_path / "native.json"
    if invalid == "bad-json":
        additions.write_text("{")
    else:
        component = dict(original)
        if invalid == "bad-hash":
            component["externalReferences"] = [
                {"type": "distribution", "hashes": [{"alg": "SHA-256", "content": "invalid"}]}
            ]
        else:
            component["version"] = "1.24"
        additions.write_text(json.dumps({"components": [component]}))
    with pytest.raises(ValueError):
        native.merge(sbom, additions)
    assert sbom.read_bytes() == before


def test_embedded_htslib_inventory_uses_loaded_versions_and_fails_broken_optional(
    tmp_path, monkeypatch
):
    native = load_script("native_sbom")
    modules = {
        "wasp2_rust": SimpleNamespace(__htslib_version__="1.19.1"),
        "pysam.libchtslib": SimpleNamespace(__file__="pysam-library.so"),
        "pysam.version": SimpleNamespace(
            __htslib_version__="incorrect-metadata",
            __samtools_version__="1.24",
            __bcftools_version__="1.24",
        ),
        "cyvcf2.cyvcf2": SimpleNamespace(__file__="cyvcf2-library.so"),
    }
    monkeypatch.setattr(native.importlib, "import_module", modules.__getitem__)
    monkeypatch.setattr(native.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(
        native,
        "_loaded_htslib_version",
        {"pysam-library.so": "1.24", "cyvcf2-library.so": "1.25"}.__getitem__,
    )
    copies = [
        component for component in native._embedded_components() if component["name"] == "htslib"
    ]
    assert {component["version"] for component in copies} == {"1.19.1", "1.24", "1.25"}
    del modules["cyvcf2.cyvcf2"]
    with pytest.raises(KeyError, match="cyvcf2"):
        native._embedded_components()


def test_loaded_htslib_calls_native_function_with_correct_signature(monkeypatch):
    native = load_script("native_sbom")

    def hts_version():
        assert hts_version.argtypes == []
        assert hts_version.restype is native.ctypes.c_char_p
        return b"1.25"

    monkeypatch.setattr(
        native.ctypes, "CDLL", lambda path: SimpleNamespace(hts_version=hts_version)
    )
    assert native._loaded_htslib_version("library.so") == "1.25"


@pytest.mark.parametrize(
    ("name", "version", "supported"),
    [
        ("htslib", "1.21.2", False),
        ("htslib", "1.21.3", True),
        ("htslib", "1.22.3", False),
        ("htslib", "1.22.4", True),
        ("htslib", "1.23.2", False),
        ("htslib", "1.23.3", True),
        ("htslib", "1.24", False),
        ("htslib", "1.24.1", True),
        ("htslib", "1.25", True),
        ("htslib", "1.19-1-g8f72310", False),
        ("htslib", "1.20.99", False),
        ("htslib", "1.25-rc1", False),
        ("htslib", "unknown", False),
        ("htslib", None, False),
        ("htslib", 1.25, False),
        ("samtools", "1.23.3", False),
        ("samtools", "1.24", False),
        ("samtools", "1.24.1", True),
        ("samtools", "1.25", True),
    ],
)
def test_native_security_baseline_boundaries(tmp_path, name, version, supported):
    native = load_script("native_sbom")
    components = [
        native._component("htslib", "1.25", "source", "manifest"),
        native._component("samtools", "1.25", "source", "manifest"),
    ]
    tested = native._component(name, "placeholder", "embedded", "native")
    tested["version"] = version
    components.append(tested)
    inventory = tmp_path / "native.json"
    inventory.write_text(json.dumps({"components": components}))
    if supported:
        native.check(inventory)
    else:
        with pytest.raises(ValueError, match="Unsupported native security baseline") as error:
            native.check(inventory)
        assert tested["bom-ref"] in str(error.value)
        assert str(version) in str(error.value)


def test_native_security_rejects_unsafe_embedded_copy_beside_current_cli(tmp_path):
    native = load_script("native_sbom")
    inventory = tmp_path / "native.json"
    components = [
        native._component("htslib", "1.25", "source", "manifest"),
        native._component("samtools", "1.25", "source", "manifest"),
        native._component("htslib", "1.25", "wasp2-rust", "hts_version()"),
        native._component("htslib", "1.24", "pysam", "hts_version()"),
        native._component("htslib", "1.19-1-g8f72310", "cyvcf2", "hts_version()"),
    ]
    inventory.write_text(json.dumps({"components": components}))
    with pytest.raises(ValueError) as error:
        native.check(inventory)
    assert "wasp2:native:pysam:htslib@1.24: 1.24" in str(error.value)
    assert "wasp2:native:cyvcf2:htslib@1.19-1-g8f72310: 1.19-1-g8f72310" in str(error.value)


@pytest.mark.parametrize("missing", ["htslib", "samtools"])
def test_native_security_requires_both_native_projects(tmp_path, missing):
    native = load_script("native_sbom")
    inventory = tmp_path / "native.json"
    inventory.write_text(
        json.dumps(
            {
                "components": [
                    native._component(name, "1.25", "source", "manifest")
                    for name in ["htslib", "samtools"]
                    if name != missing
                ]
            }
        )
    )
    with pytest.raises(ValueError, match="missing required HTSlib or Samtools"):
        native.check(inventory)
