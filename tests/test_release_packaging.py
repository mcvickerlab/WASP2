"""Publication guards must fail before changing files or uploading artifacts."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path
from urllib.error import HTTPError

import pytest

ROOT = Path(__file__).resolve().parents[1]


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
