"""Small deterministic mapping fixtures; no external benchmark checkout required."""

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import pysam


@dataclass(frozen=True)
class Variant:
    start: int
    ref: str
    alt: str


def quickbench_snv_variants():
    return [Variant(110, "A", "C"), Variant(210, "A", "G"), Variant(220, "A", "T")]


def quickbench_indel_variants():
    return [Variant(110, "A", "ATT")]


def write_bed(path: Path, variants):
    path.write_text(
        "".join(
            f"chr1\t{v.start}\t{v.start + len(v.ref)}\t{v.ref}\t{v.alt}\t{v.ref}|{v.alt}\n"
            for v in variants
        )
    )


def write_pair(path: Path, name: str, reverse_mate: bool):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for mate, position in enumerate((100, 200)):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.flag = (99 if mate == 0 else 147) if reverse_mate else (67 if mate == 0 else 131)
            read.reference_id = read.next_reference_id = 0
            read.reference_start = position
            read.next_reference_start = 200 if mate == 0 else 100
            read.template_length = 150 if mate == 0 else -150
            read.mapping_quality = 60
            read.cigarstring = "50M"
            read.query_sequence = "A" * 50
            read.query_qualities = list(range(50))
            bam.write(read)
    pysam.index(str(path))


def write_synthetic_bam(path):
    write_pair(path, "pairR", True)


def write_synthetic_bam_indel(path):
    write_pair(path, "pairI", False)


class CanonicalFastqRecord(NamedTuple):
    name: str
    mate: int
    sequence: str
    quality: str


def iter_fastq(path):
    with pysam.FastxFile(str(path)) as records:
        for record in records:
            yield record.name, record.sequence, record.quality


def fastq_counter(*paths):
    result = Counter()
    for mate, path in enumerate(paths, start=1):
        for name, sequence, quality in iter_fastq(path):
            result[CanonicalFastqRecord(name.split("_WASP_", 1)[0], mate, sequence, quality)] += 1
    return result


def counter_diff(left, right):
    return list((left - right).elements()), list((right - left).elements())
