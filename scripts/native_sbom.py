#!/usr/bin/env python3
"""Record native source identities and independently loaded HTSlib copies."""

import argparse
import ctypes
import importlib
import importlib.util
import json
import os
import re
import tempfile
from pathlib import Path
from urllib.parse import quote, urlsplit


def _version(value: object) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError("native component version must be a nonempty string")
    return value


def _component(name: str, version: object, origin: str, evidence: str) -> dict:
    version = _version(version)
    return {
        "type": "library" if name == "htslib" else "application",
        "name": name,
        "version": version,
        "bom-ref": f"wasp2:native:{origin}:{name}@{quote(version, safe='')}",
        "purl": f"pkg:generic/{name}@{quote(version, safe='')}",
        "properties": [
            {"name": "wasp2:native-origin", "value": origin},
            {"name": "wasp2:version-evidence", "value": evidence},
        ],
    }


def _source_components(manifest: Path) -> list[dict]:
    document = json.loads(manifest.read_text())
    records = document.get("components") if isinstance(document, dict) else None
    if not isinstance(records, list):
        raise ValueError("native manifest requires a components list")
    components = []
    names = set()
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("native manifest component must be an object")
        name = record.get("name")
        if (
            not isinstance(name, str)
            or name not in {"htslib", "samtools", "bcftools"}
            or name in names
        ):
            raise ValueError("native manifest has an unexpected or duplicate component")
        names.add(name)
        digest = record.get("sha256")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-fA-F]{64}", digest) is None:
            raise ValueError(f"{name}: invalid source SHA-256")
        url = record.get("url")
        if not isinstance(url, str) or urlsplit(url).scheme != "https" or not urlsplit(url).netloc:
            raise ValueError(f"{name}: source URL must be HTTPS")
        component = _component(name, record.get("version"), "source", "source manifest")
        component["externalReferences"] = [
            {
                "type": "distribution",
                "url": url,
                "hashes": [{"alg": "SHA-256", "content": digest.lower()}],
            }
        ]
        components.append(component)
    if names != {"htslib", "samtools", "bcftools"}:
        raise ValueError("native manifest must inventory htslib, samtools, and bcftools")
    return components


def _loaded_htslib_version(path: str) -> str:
    library = ctypes.CDLL(path)
    function = library.hts_version
    function.argtypes = []
    function.restype = ctypes.c_char_p
    value = function()
    if not value:
        raise ValueError(f"{path}: hts_version returned no version")
    return _version(value.decode("utf-8"))


def _embedded_components() -> list[dict]:
    rust = importlib.import_module("wasp2_rust")
    pysam_library = importlib.import_module("pysam.libchtslib")
    pysam_version = importlib.import_module("pysam.version")
    components = [
        _component("htslib", rust.__htslib_version__, "wasp2-rust", "hts_version()"),
        _component(
            "htslib",
            _loaded_htslib_version(pysam_library.__file__),
            "pysam",
            "loaded pysam.libchtslib hts_version()",
        ),
        _component("samtools", pysam_version.__samtools_version__, "pysam", "pysam.version"),
        _component("bcftools", pysam_version.__bcftools_version__, "pysam", "pysam.version"),
    ]
    if importlib.util.find_spec("cyvcf2") is not None:
        cyvcf2_library = importlib.import_module("cyvcf2.cyvcf2")
        components.append(
            _component(
                "htslib",
                _loaded_htslib_version(cyvcf2_library.__file__),
                "cyvcf2",
                "loaded cyvcf2.cyvcf2 hts_version()",
            )
        )
    return components


def collect(manifest: Path) -> dict:
    components = _source_components(manifest)
    components.extend(_embedded_components())
    return {"components": components}


def check(native: Path) -> None:
    """Enforce supported native security baselines independently of image scans."""
    # HTSlib security fixes and maintained branches:
    # https://github.com/samtools/htslib/releases/tag/1.25
    # Samtools minimum supported baseline is 1.24.1; this does not assert that
    # every earlier release is vulnerable (older branches have backports):
    # https://github.com/samtools/samtools/security/advisories/GHSA-vh3p-9c59-j3hv
    document = json.loads(native.read_text())
    components = document.get("components") if isinstance(document, dict) else None
    if not isinstance(components, list) or not components:
        raise ValueError("native security check requires components")
    failures = []
    seen = set()
    for component in components:
        if not isinstance(component, dict):
            raise ValueError("native component must be an object")
        name = component.get("name")
        if name not in {"htslib", "samtools"}:
            continue
        seen.add(name)
        version = component.get("version")
        reference = component.get("bom-ref")
        # Accept official release versions and git-describe suffixes, never
        # infer a patched release from an unknown or prerelease version string.
        match = (
            re.fullmatch(r"(\d+)\.(\d+)(?:\.(\d+))?(?:-\d+-g[0-9a-fA-F]+)?", version)
            if isinstance(version, str)
            else None
        )
        supported = False
        if match is not None:
            major, minor, patch = (int(value or 0) for value in match.groups())
            if major == 1:
                if name == "htslib":
                    minimum_patch = {21: 3, 22: 4, 23: 3, 24: 1}.get(minor)
                    supported = minor >= 25 or (
                        minimum_patch is not None and patch >= minimum_patch
                    )
                else:
                    supported = minor >= 25 or (minor == 24 and patch >= 1)
        if not supported or not isinstance(reference, str) or not reference:
            failures.append(f"{reference or name}: {version}")
    if seen != {"htslib", "samtools"}:
        failures.append("missing required HTSlib or Samtools inventory")
    if failures:
        raise ValueError("Unsupported native security baseline:\n" + "\n".join(failures))


def merge(sbom: Path, native: Path) -> None:
    inventory = json.loads(sbom.read_text())
    additions = json.loads(native.read_text())
    if not isinstance(inventory, dict) or inventory.get("bomFormat") != "CycloneDX":
        raise ValueError("expected a CycloneDX SBOM")
    existing = inventory.get("components")
    records = additions.get("components") if isinstance(additions, dict) else None
    if not isinstance(existing, list) or not isinstance(records, list) or not records:
        raise ValueError("SBOM and native inventory require nonempty components lists")
    for component in records:
        if not isinstance(component, dict) or component.get("type") not in {
            "library",
            "application",
        }:
            raise ValueError("invalid native component")
        _version(component.get("version"))
        if not isinstance(component.get("name"), str) or not component["name"]:
            raise ValueError("native component is missing name")
        for reference in component.get("externalReferences", []):
            for digest in reference.get("hashes", []):
                if (
                    digest.get("alg") != "SHA-256"
                    or re.fullmatch(r"[0-9a-fA-F]{64}", digest.get("content", "")) is None
                ):
                    raise ValueError("invalid native source SHA-256")
    by_reference = {}
    for component in existing + records:
        if not isinstance(component, dict):
            raise ValueError("SBOM component must be an object")
        reference = component.get("bom-ref")
        if not isinstance(reference, str) or not reference:
            raise ValueError("SBOM component is missing bom-ref")
        if reference in by_reference and component != by_reference[reference]:
            raise ValueError(f"conflicting component identity: {reference}")
        by_reference[reference] = component
    inventory["components"] = list(by_reference.values())
    # Replace only after every input and component identity has been validated.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=sbom.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(inventory, handle, indent=2)
            handle.write("\n")
        os.chmod(temporary, sbom.stat().st_mode & 0o777)
        os.replace(temporary, sbom)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    collector = commands.add_parser("collect")
    collector.add_argument("manifest", type=Path)
    checker = commands.add_parser("check")
    checker.add_argument("native", type=Path)
    merger = commands.add_parser("merge")
    merger.add_argument("sbom", type=Path)
    merger.add_argument("native", type=Path)
    args = parser.parse_args()
    if args.command == "collect":
        print(json.dumps(collect(args.manifest), indent=2))
    elif args.command == "check":
        check(args.native)
    else:
        merge(args.sbom, args.native)


if __name__ == "__main__":
    main()
