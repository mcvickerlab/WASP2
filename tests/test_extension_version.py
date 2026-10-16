import importlib.metadata

import pytest

import wasp2_rust


def test_rust_extension_version_matches_distribution() -> None:
    assert wasp2_rust.__version__ == importlib.metadata.version("wasp2")


@pytest.fixture
def release_root(tmp_path, monkeypatch):
    """Exercise release synchronization without changing checkout metadata."""
    from scripts import release_version

    paths = [
        "rust/Cargo.toml",
        "rust/Cargo.lock",
        "pixi.toml",
        "CITATION.cff",
        "Dockerfile",
        "Singularity.def",
        "README.md",
        "docs/CONTAINER_USAGE.md",
        "docs/source/installation.rst",
        "docs/source/conf.py",
        "src/wasp2/__init__.py",
        "galaxy/tools/wasp2/macros.xml",
    ]
    for relative in paths:
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(release_version.read(relative), encoding="utf-8")
    monkeypatch.setattr(release_version, "ROOT", tmp_path)
    return release_version


@pytest.mark.parametrize("root_version", ["current", "missing", "9.8.7"])
def test_release_version_checks_root_lock_package(release_root, root_version):
    import re

    version = release_root.cargo_version()
    lock = release_root.read("rust/Cargo.lock")
    # Leave dependency versions intact: only the root package should determine
    # the lockfile's release-version compatibility.
    pattern = rf'(?m)^name = "wasp2"\nversion = "{re.escape(version)}"'
    replacement = (
        'name = "different-package"\nversion = "9.8.7"'
        if root_version == "missing"
        else f'name = "wasp2"\nversion = "{version if root_version == "current" else root_version}"'
    )
    updated, count = re.subn(pattern, lambda _match: replacement, lock)
    assert count == 1
    release_root.write("rust/Cargo.lock", updated)
    errors = release_root.check_version()
    if root_version == "current":
        assert errors == []
    else:
        assert errors == ["rust/Cargo.lock root package version differs from rust/Cargo.toml"]


@pytest.mark.parametrize("existing_date", [None, "2020-01-01"])
def test_setting_release_version_inserts_or_updates_citation_date(release_root, existing_date):
    import re

    citation = re.sub(r"(?m)^date-released:.*\n?", "", release_root.read("CITATION.cff"))
    if existing_date is not None:
        citation += f'\ndate-released: "{existing_date}"\n'
    release_root.write("CITATION.cff", citation)
    release_root.set_version("9.8.7", "2027-02-03")
    updated = release_root.read("CITATION.cff")
    assert re.findall(r"(?m)^date-released:.*$", updated) == ['date-released: "2027-02-03"']
    assert 'version: "9.8.7"' in updated
    assert release_root.cargo_version() == "9.8.7"
    assert release_root.check_version() == []
