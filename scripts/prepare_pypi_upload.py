#!/usr/bin/env python3
"""Stage missing PyPI artifacts, refusing to overwrite different published bytes."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.error
import urllib.request
from pathlib import Path


def stage(packages: Path, destination: Path, version: str) -> int:
    try:
        with urllib.request.urlopen(
            f"https://pypi.org/pypi/wasp2/{version}/json", timeout=30
        ) as response:
            published = json.load(response)["urls"]
    except urllib.error.HTTPError as error:
        if error.code != 404:
            raise
        published = []
    hashes = {item["filename"]: item["digests"]["sha256"] for item in published}
    artifacts = sorted(p for p in packages.iterdir() if p.name.endswith((".whl", ".tar.gz")))
    if not artifacts:
        raise ValueError("No package artifacts supplied")
    missing = []
    for artifact in artifacts:
        expected = hashes.get(artifact.name)
        if expected is None:
            missing.append(artifact)
        elif hashlib.sha256(artifact.read_bytes()).hexdigest() != expected:
            raise ValueError(
                f"Published artifact differs: {artifact.name}; do not overwrite a release"
            )
    # Check all hashes before creating any upload payload.
    destination.mkdir(parents=True, exist_ok=False)
    for artifact in missing:
        shutil.copy2(artifact, destination / artifact.name)
    return len(missing)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packages", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("version")
    args = parser.parse_args()
    count = stage(args.packages, args.destination, args.version)
    print(f"Staged {count} missing artifacts")
    # GitHub Actions reads this output; stdout remains useful locally.
    import os

    if output := os.environ.get("GITHUB_OUTPUT"):
        with open(output, "a") as handle:
            handle.write(f"count={count}\n")


if __name__ == "__main__":
    main()
