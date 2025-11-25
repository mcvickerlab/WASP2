"""Helpers to build synthetic mapping test data."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Tuple

import pysam


def write_bams_from_real(
    real_bam: Path,
    out_dir: Path,
    n_pairs: int = 20000,
    moved_fraction: float = 0.2,
) -> Tuple[Path, Path]:
    """Create synthetic original/remapped BAMs for testing mapping filter."""
    out_dir.mkdir(parents=True, exist_ok=True)
    orig_path = out_dir / "orig.bam"
    remap_path = out_dir / "remap.bam"

    with pysam.AlignmentFile(real_bam, "rb") as bam_in:
        header = bam_in.header.to_dict()
        with pysam.AlignmentFile(orig_path, "wb", header=header) as bam_orig, pysam.AlignmentFile(
            remap_path, "wb", header=header
        ) as bam_remap:
            pairs_written = 0
            for read in bam_in.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                # Take only first of each pair
                if read.is_read2:
                    continue
                try:
                    mate = bam_in.mate(read) if read.is_paired else None
                except ValueError:
                    # mate not found in this small test BAM; skip
                    continue
                if mate and mate.is_unmapped:
                    continue

                # assign synthetic WASP name encoding
                orig_pos = (read.reference_start, mate.reference_start if mate else -1)
                total_copies = 2 if mate else 1
                read.query_name = f"{read.query_name}_WASP_{orig_pos[0]}_{orig_pos[1]}_{total_copies}"
                if mate:
                    mate.query_name = read.query_name

                bam_orig.write(read)
                if mate:
                    bam_orig.write(mate)

                # clone remapped; optionally move position
                def moved(r):
                    return random.random() < moved_fraction

                r1 = read.to_string().split("\t")
                if moved(read):
                    # shift by +1000 for synthetic move
                    r1[3] = str(read.reference_start + 1000 + 1)  # 1-based POS in SAM
                bam_remap.write(pysam.AlignedSegment.fromstring("\t".join(r1), bam_remap.header))

                if mate:
                    r2 = mate.to_string().split("\t")
                    if moved(mate):
                        r2[3] = str(mate.reference_start + 1000 + 1)
                    bam_remap.write(pysam.AlignedSegment.fromstring("\t".join(r2), bam_remap.header))

                pairs_written += 1
                if pairs_written >= n_pairs:
                    break

    return orig_path, remap_path


def collect_qnames(bam_path: Path) -> set[str]:
    names: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            names.add(read.query_name)
    return names
