#!/usr/bin/env bash
# Build the command-line HTS tools independently of distro HTSlib packages.
set -euo pipefail

prefix=/opt/wasp2/hts-tools
version=1.25
work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT
cd "$work_dir"

names=(htslib samtools bcftools)
checksums=(
    158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730
    8bbd2b3cb447953149a74ab0b00e3ec4e88a789a27b8373a31fa6b5b2dee02b1
    f96e583d2315998645bb56756f1c8f230c489512db2655b0f623198ebc8a2d99
)

for i in "${!names[@]}"; do
    name=${names[$i]}
    archive="$name-$version.tar.bz2"
    url="https://github.com/samtools/$name/releases/download/$version/$archive"
    curl --fail --location --retry 5 --retry-all-errors --output "$archive" "$url"
    printf '%s  %s\n' "${checksums[$i]}" "$archive" | sha256sum --check --strict
    tar -xjf "$archive"
done

jobs=$(getconf _NPROCESSORS_ONLN)
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="$prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

(
    cd "htslib-$version"
    ./configure --prefix="$prefix" --enable-libcurl --enable-s3 --enable-gcs \
        --with-libdeflate
    make -j "$jobs"
    make install
)
(
    cd "samtools-$version"
    ./configure --prefix="$prefix" --with-htslib="$prefix" --with-ncursesw
    make -j "$jobs"
    make install
)
(
    cd "bcftools-$version"
    ./configure --prefix="$prefix" --with-htslib="$prefix" \
        --enable-bcftools-plugins --enable-libgsl --with-cblas=gslcblas
    make -j "$jobs"
    make install
)

for name in "${names[@]}"; do
    mkdir -p "$prefix/share/licenses/$name"
    install -m 644 "$name-$version/LICENSE" "$prefix/share/licenses/$name/LICENSE"
done

# Keep the exact source identities available to the runtime SBOM generator.
python3 - "$version" "${checksums[@]}" <<'PY'
import json
import sys
from pathlib import Path

version, *checksums = sys.argv[1:]
components = [
    {
        "name": name,
        "version": version,
        "sha256": checksum,
        "url": (
            f"https://github.com/samtools/{name}/releases/download/"
            f"{version}/{name}-{version}.tar.bz2"
        ),
    }
    for name, checksum in zip(("htslib", "samtools", "bcftools"), checksums, strict=True)
]
Path("/opt/wasp2/native-sources.json").write_text(
    json.dumps({"components": components}, indent=2) + "\n"
)
PY

"$prefix/bin/samtools" --version
"$prefix/bin/bcftools" --version
"$prefix/bin/tabix" --version
"$prefix/bin/bcftools" plugin --list | grep -Fx fill-tags
