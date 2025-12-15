from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Counter as CounterType, Iterable, Iterator, Tuple


@dataclass(frozen=True)
class CanonicalFastqRecord:
    orig_name: str
    mate: int
    sequence: str
    qualities: str


def _open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def iter_fastq(path: Path) -> Iterator[Tuple[str, str, str]]:
    """Yield (name, seq, qual) from a FASTQ (optionally gzipped)."""
    with _open_maybe_gzip(path) as f:
        while True:
            name = f.readline()
            if not name:
                return
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                return
            if not name.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield name[1:].strip().split()[0], seq.strip(), qual.strip()


def canonicalize_wasp_name(name: str) -> Tuple[str, int]:
    """Return (orig_name, mate) from a WASP FASTQ qname.

    Accepts names like:
    - readX_WASP_100_200_1_2/1
    - readX_WASP_100_200_1_1_3_4/2
    """
    if name.endswith("/1"):
        mate = 1
        core = name[:-2]
    elif name.endswith("/2"):
        mate = 2
        core = name[:-2]
    else:
        raise ValueError(f"FASTQ name missing /1 or /2 mate suffix: {name}")

    split = core.split("_WASP_", 1)
    orig_name = split[0]
    return orig_name, mate


def fastq_counter(r1_path: Path, r2_path: Path) -> CounterType[CanonicalFastqRecord]:
    c: CounterType[CanonicalFastqRecord] = Counter()
    for path in (r1_path, r2_path):
        for name, seq, qual in iter_fastq(path):
            orig, mate = canonicalize_wasp_name(name)
            c[CanonicalFastqRecord(orig, mate, seq, qual)] += 1
    return c


def counter_diff(
    a: CounterType[CanonicalFastqRecord],
    b: CounterType[CanonicalFastqRecord],
) -> tuple[list[tuple[CanonicalFastqRecord, int]], list[tuple[CanonicalFastqRecord, int]]]:
    """Return (only_in_a, only_in_b) with multiplicities."""
    only_in_a: list[tuple[CanonicalFastqRecord, int]] = []
    only_in_b: list[tuple[CanonicalFastqRecord, int]] = []

    keys = set(a) | set(b)
    for k in sorted(keys, key=lambda r: (r.orig_name, r.mate, r.sequence, r.qualities)):
        da = a.get(k, 0) - b.get(k, 0)
        if da > 0:
            only_in_a.append((k, da))
        db = b.get(k, 0) - a.get(k, 0)
        if db > 0:
            only_in_b.append((k, db))
    return only_in_a, only_in_b

