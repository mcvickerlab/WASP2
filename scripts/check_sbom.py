#!/usr/bin/env python3
"""Check that a CycloneDX SBOM actually inventories required components."""

import argparse
import json
from pathlib import Path


def check(path: Path, required: list[str]) -> None:
    inventory = json.loads(path.read_text())
    names = {item["name"].lower().replace("_", "-") for item in inventory.get("components", [])}
    missing = {name.lower().replace("_", "-") for name in required} - names
    if missing:
        raise ValueError(f"{path}: missing SBOM components: {', '.join(sorted(missing))}")
    print(f"Verified {len(names)} inventoried components in {path.name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("required", nargs="+")
    args = parser.parse_args()
    check(args.path, args.required)
