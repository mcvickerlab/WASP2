#!/usr/bin/env python3
"""Synchronize and validate the active WASP2 release version."""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEMVER = r"[0-9]+\.[0-9]+\.[0-9]+(?:[a-z0-9.-]+)?"


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, content: str) -> None:
    (ROOT / path).write_text(content, encoding="utf-8")


def toml_section_version(path: str, section: str) -> str:
    content = read(path)
    section_match = re.search(
        rf"^\[{re.escape(section)}\]\s*$([\s\S]*?)(?=^\[|\Z)", content, flags=re.MULTILINE
    )
    if section_match is None:
        raise RuntimeError(f"missing [{section}] section in {path}")
    version_match = re.search(
        rf'^version\s*=\s*"({SEMVER})"', section_match.group(1), flags=re.MULTILINE
    )
    if version_match is None:
        raise RuntimeError(f"missing version in [{section}] section of {path}")
    return version_match.group(1)


def cargo_version() -> str:
    return toml_section_version("rust/Cargo.toml", "package")


def replace_one(path: str, pattern: str, replacement: str) -> None:
    content = read(path)
    updated, count = re.subn(pattern, replacement, content, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"expected one version field in {path}, found {count}")
    write(path, updated)


def release_pin_files() -> list[Path]:
    roots = [ROOT / "pipelines", ROOT / "galaxy"]
    suffixes = {".nf", ".config", ".yml", ".yaml", ".xml"}
    return sorted(
        path
        for root in roots
        for path in root.rglob("*")
        if path.is_file() and path.suffix in suffixes
    )


def set_version(version: str, release_date: str) -> None:
    if re.fullmatch(SEMVER, version) is None:
        raise RuntimeError(f"invalid release version: {version}")

    replace_one("rust/Cargo.toml", rf'^(version\s*=\s*)"{SEMVER}"', rf'\1"{version}"')
    replace_one(
        "rust/Cargo.lock",
        rf'(^name\s*=\s*"wasp2"\s*\nversion\s*=\s*)"{SEMVER}"',
        rf'\1"{version}"',
    )
    replace_one("pixi.toml", rf'^(version\s*=\s*)"{SEMVER}"', rf'\1"{version}"')
    replace_one("CITATION.cff", rf'^(version:\s*)"{SEMVER}"', rf'\1"{version}"')
    replace_one(
        "CITATION.cff",
        r'^(date-released:\s*)"[0-9]{4}-[0-9]{2}-[0-9]{2}"',
        rf'\1"{release_date}"',
    )
    replace_one("Dockerfile", rf"^(ARG VERSION=){SEMVER}", rf"\g<1>{version}")
    replace_one(
        "Singularity.def",
        rf"^(From:\s*ghcr\.io/mcvickerlab/wasp2:){SEMVER}",
        rf"\g<1>{version}",
    )
    replace_one("Singularity.def", rf"^(\s*Version\s+){SEMVER}", rf"\g<1>{version}")

    direct_files = ["README.md", "docs/CONTAINER_USAGE.md", "docs/source/installation.rst"]
    for path in direct_files:
        content = re.sub(rf"(ghcr\.io/mcvickerlab/wasp2:){SEMVER}", rf"\g<1>{version}", read(path))
        write(path, content)

    patterns = [
        (rf"(ghcr\.io/mcvickerlab/wasp2:){SEMVER}", rf"\g<1>{version}"),
        (rf"(wasp2==){SEMVER}", rf"\g<1>{version}"),
        (rf"(bioconda::wasp2=){SEMVER}", rf"\g<1>{version}"),
        (rf"(^\s*wasp2:\s*){SEMVER}(\s*$)", rf"\g<1>{version}\2"),
        (rf"(wasp2_container_version\s*=\s*'){SEMVER}(')", rf"\g<1>{version}\2"),
        (rf'(<token name="@TOOL_VERSION@">){SEMVER}(</token>)', rf"\g<1>{version}\2"),
    ]
    for path in release_pin_files():
        content = path.read_text(encoding="utf-8")
        updated = content
        for pattern, replacement in patterns:
            updated = re.sub(pattern, replacement, updated, flags=re.MULTILINE)
        if updated != content:
            path.write_text(updated, encoding="utf-8")


def check_version() -> list[str]:
    expected = cargo_version()
    errors: list[str] = []

    pixi_version = toml_section_version("pixi.toml", "project")
    if pixi_version != expected:
        errors.append(f"pixi.toml: {pixi_version} != {expected}")

    checks = {
        "CITATION.cff": rf'^version:\s*"({SEMVER})"',
        "Dockerfile": rf"^ARG VERSION=({SEMVER})",
        "Singularity.def (image)": rf"^From:\s*ghcr\.io/mcvickerlab/wasp2:({SEMVER})",
        "Singularity.def (label)": rf"^\s*Version\s+({SEMVER})",
        "galaxy/tools/wasp2/macros.xml": rf'<token name="@TOOL_VERSION@">({SEMVER})</token>',
    }
    for label, pattern in checks.items():
        path = label.split(" (")[0]
        match = re.search(pattern, read(path), flags=re.MULTILINE)
        found = match.group(1) if match else "missing"
        if found != expected:
            errors.append(f"{label}: {found} != {expected}")

    active_patterns = [
        rf"ghcr\.io/mcvickerlab/wasp2:({SEMVER})",
        rf"wasp2==({SEMVER})",
        rf"bioconda::wasp2=({SEMVER})",
        rf"^\s*wasp2:\s*({SEMVER})\s*$",
        rf"wasp2_container_version\s*=\s*'({SEMVER})'",
    ]
    direct_files = [
        ROOT / "README.md",
        ROOT / "docs/CONTAINER_USAGE.md",
        ROOT / "docs/source/installation.rst",
    ]
    for path in [*direct_files, *release_pin_files()]:
        content = path.read_text(encoding="utf-8")
        for pattern in active_patterns:
            for match in re.finditer(pattern, content, flags=re.MULTILINE):
                if match.group(1) != expected:
                    errors.append(f"{path.relative_to(ROOT)}: {match.group(1)} != {expected}")

    init_content = read("src/wasp2/__init__.py")
    if 'version("wasp2")' not in init_content:
        errors.append("src/wasp2/__init__.py does not derive its version from package metadata")
    if "tomllib.load" not in read("docs/source/conf.py"):
        errors.append("docs/source/conf.py does not derive its version from Cargo.toml")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true", help="validate all active version pins")
    group.add_argument("--set", metavar="VERSION", help="update all active version pins")
    parser.add_argument("--date", default=dt.date.today().isoformat(), help="release date")
    args = parser.parse_args()

    if args.set:
        set_version(args.set, args.date)

    errors = check_version()
    if errors:
        print("Version consistency failed:", file=sys.stderr)
        for error in sorted(set(errors)):
            print(f"  - {error}", file=sys.stderr)
        return 1

    print(f"All active release surfaces use WASP2 {cargo_version()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
