from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

import pysam


@dataclass(frozen=True)
class Variant:
    chrom: str
    start: int  # 0-based
    ref: str
    alt: str
    genotype: str  # e.g. "A|G" or "0|1"

    @property
    def stop(self) -> int:
        return self.start + len(self.ref)

    def to_bed_row(self) -> str:
        return f"{self.chrom}\t{self.start}\t{self.stop}\t{self.ref}\t{self.alt}\t{self.genotype}\n"


def write_bed(path: Path, variants: Sequence[Variant]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for v in variants:
            f.write(v.to_bed_row())


def _make_seq(length: int, replacements: dict[int, str]) -> str:
    seq = ["A"] * length
    for idx, base in replacements.items():
        seq[idx] = base
    return "".join(seq)


def _write_pair(
    out: pysam.AlignmentFile,
    *,
    header: pysam.AlignmentHeader,
    qname: str,
    r1_start: int,
    r2_start: int,
    r1_seq: str,
    r2_seq: str,
    mapq: int = 60,
) -> None:
    # Flags chosen to keep everything simple: paired, proper_pair, forward strand.
    r1_flag = 0x1 | 0x2 | 0x40
    r2_flag = 0x1 | 0x2 | 0x80

    r1 = pysam.AlignedSegment(header)
    r1.query_name = qname
    r1.flag = r1_flag
    r1.reference_id = 0
    r1.reference_start = r1_start
    r1.mapping_quality = mapq
    r1.cigarstring = f"{len(r1_seq)}M"
    r1.next_reference_id = 0
    r1.next_reference_start = r2_start
    r1.template_length = (r2_start + len(r2_seq)) - r1_start
    r1.query_sequence = r1_seq
    r1.query_qualities = pysam.qualitystring_to_array("I" * len(r1_seq))

    r2 = pysam.AlignedSegment(header)
    r2.query_name = qname
    r2.flag = r2_flag
    r2.reference_id = 0
    r2.reference_start = r2_start
    r2.mapping_quality = mapq
    r2.cigarstring = f"{len(r2_seq)}M"
    r2.next_reference_id = 0
    r2.next_reference_start = r1_start
    r2.template_length = -r1.template_length
    r2.query_sequence = r2_seq
    r2.query_qualities = pysam.qualitystring_to_array("I" * len(r2_seq))

    out.write(r1)
    out.write(r2)


def write_synthetic_bam(path: Path) -> None:
    """Write a tiny coordinate-sorted paired-end BAM for quickbench.

    This dataset is designed to exercise:
    - A pair where the read's alleles DON'T match either haplotype (→ 2 remap copies)
    - A pair where one haplotype matches the original (→ 1 remap copy)
    - A pair with no overlaps (→ kept, no remap reads)
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.0", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
        }
    )

    unsorted_path = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as out:
        # Variants used by quickbench_snv_variants():
        # - chr1:100 A/G (hap1=A, hap2=G)
        # - chr1:110 C/T (hap1=C, hap2=T)
        #
        # Pair A: R1 overlaps both variants and is a "mix" (A at 100, T at 110),
        # which doesn't match hap1 (A,C) or hap2 (G,T) → 2 remap copies.
        r1_a = _make_seq(50, {10: "A", 20: "T"})  # start 90 → offsets 10 and 20
        r2_a = "G" * 50
        _write_pair(
            out,
            header=header,
            qname="pairA",
            r1_start=90,
            r2_start=400,
            r1_seq=r1_a,
            r2_seq=r2_a,
        )

        # Pair B: R1 overlaps both variants and matches hap1 (A at 100, C at 110)
        # → only the hap2 remap copy differs.
        r1_b = _make_seq(50, {10: "A", 20: "C"})  # start 91 → offsets 9 and 19 would shift
        r2_b = "T" * 50
        _write_pair(
            out,
            header=header,
            qname="pairB",
            r1_start=90,
            r2_start=401,
            r1_seq=r1_b,
            r2_seq=r2_b,
        )

        # Pair C: no overlap with variants (kept, no remap reads).
        r1_c = "C" * 50
        r2_c = "A" * 50
        _write_pair(
            out,
            header=header,
            qname="pairC",
            r1_start=700,
            r2_start=760,
            r1_seq=r1_c,
            r2_seq=r2_c,
        )

    # Coordinate-sort + index for downstream convenience.
    pysam.sort("-o", str(path), str(unsorted_path))
    pysam.index(str(path))
    unsorted_path.unlink(missing_ok=True)


def quickbench_snv_variants() -> List[Variant]:
    return [
        Variant(chrom="chr1", start=100, ref="A", alt="G", genotype="A|G"),
        Variant(chrom="chr1", start=110, ref="C", alt="T", genotype="C|T"),
    ]


def quickbench_indel_variants() -> List[Variant]:
    # Simple 2bp insertion represented with an anchor base (VCF-style):
    # ref="A", alt="ACG" → +2 bp net insertion.
    return [Variant(chrom="chr1", start=100, ref="A", alt="ACG", genotype="A|ACG")]


def write_synthetic_bam_indel(path: Path) -> None:
    """Write a tiny coordinate-sorted paired-end BAM for INDEL quickbench.

    Includes:
    - One pair overlapping a +2bp insertion (hap2 differs; hap1 matches original)
    - One pair with no overlaps
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.0", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
        }
    )

    unsorted_path = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as out:
        # Pair I: R1 overlaps insertion at ref pos 100 (read starts at 90).
        # Original sequence matches hap1 ("A" at the anchor base), so only hap2 differs.
        r1_i = ("T" * 5) + ("A" * 40) + ("C" * 5)  # length 50, distinct ends
        r2_i = "G" * 50
        _write_pair(
            out,
            header=header,
            qname="pairI",
            r1_start=90,
            r2_start=400,
            r1_seq=r1_i,
            r2_seq=r2_i,
        )

        # Pair N: no overlap with the variant.
        r1_n = "C" * 50
        r2_n = "A" * 50
        _write_pair(
            out,
            header=header,
            qname="pairN",
            r1_start=700,
            r2_start=760,
            r1_seq=r1_n,
            r2_seq=r2_n,
        )

    pysam.sort("-o", str(path), str(unsorted_path))
    pysam.index(str(path))
    unsorted_path.unlink(missing_ok=True)
