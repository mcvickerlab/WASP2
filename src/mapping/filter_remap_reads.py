import os
import subprocess
import timeit
from typing import Optional

import pysam

try:
    from wasp2_rust import filter_bam_wasp
except ImportError:
    filter_bam_wasp = None

def _filter_bam_wasp_python(
    to_remap_bam: str,
    remapped_bam: str,
    filt_out_bam: str,
    keep_read_file: Optional[str] = None,
    threads: int = 1,
    same_locus_slop: int = 0,
) -> None:
    """Python fallback for WASP filtering (slower than Rust but works)."""
    import timeit
    start_time = timeit.default_timer()

    # Track expected positions and remaining remapped copies
    keep_set = set()
    pos_map = {}  # read_name -> (pos1, pos2)
    remaining = {}  # read_name -> count remaining
    removed_moved = 0

    # Phase 1: Scan remapped BAM to identify reads that pass WASP filter
    print(f"Phase 1: Scanning remapped BAM...")
    with pysam.AlignmentFile(remapped_bam, "rb", threads=threads) as bam:
        for rec in bam:
            if rec.is_unmapped or not rec.is_proper_pair or rec.is_secondary or rec.is_supplementary:
                continue

            # Parse _WASP_ encoded name: orig_name_WASP_pos1_pos2_copynum_total
            qname = rec.query_name
            if "_WASP_" not in qname:
                continue

            parts = qname.split("_WASP_")
            if len(parts) != 2:
                continue

            orig_name = parts[0]
            suffix_parts = parts[1].split("_")
            if len(suffix_parts) < 4:
                continue

            try:
                pos1 = int(suffix_parts[0])
                pos2 = int(suffix_parts[1])
                total = int(suffix_parts[3])
            except ValueError:
                continue

            # Initialize tracking for this read
            if orig_name not in pos_map:
                pos_map[orig_name] = (pos1, pos2)
                remaining[orig_name] = total
                keep_set.add(orig_name)
            elif orig_name not in keep_set:
                continue

            # Decrement remaining count
            if orig_name in remaining:
                remaining[orig_name] -= 1

            # Check if remapped position matches original (mate order agnostic)
            # NOTE: WASP encoding uses 0-based coords, rec.pos() returns 0-based position (same as reference_start)
            # BUT: The Rust code uses rec.pos() which in rust-htslib is 0-based
            # In pysam, rec.reference_start is 0-based, so we should compare directly
            rec_pos = rec.reference_start  # 0-based
            mate_pos = rec.next_reference_start  # 0-based
            expect_pos, expect_mate = pos_map[orig_name]  # Already 0-based from WASP encoding

            if same_locus_slop == 0:
                # Strict matching for SNPs
                matches = (
                    (rec_pos == expect_pos and mate_pos == expect_mate)
                    or (rec_pos == expect_mate and mate_pos == expect_pos)
                )
            else:
                # Allow slop tolerance for indels
                pos_diff1 = abs(rec_pos - expect_pos)
                mate_diff1 = abs(mate_pos - expect_mate)
                pos_diff2 = abs(rec_pos - expect_mate)
                mate_diff2 = abs(mate_pos - expect_pos)

                matches = (
                    (pos_diff1 <= same_locus_slop and mate_diff1 <= same_locus_slop)
                    or (pos_diff2 <= same_locus_slop and mate_diff2 <= same_locus_slop)
                )

            if not matches:
                keep_set.discard(orig_name)
                remaining.pop(orig_name, None)
                removed_moved += 1

    # Remove reads with missing counts
    missing_count = len(remaining)
    if missing_count > 0:
        for name in list(remaining.keys()):
            keep_set.discard(name)
        removed_moved += missing_count

    print(f"Kept {len(keep_set)} reads, removed {removed_moved} reads (Phase 1: {timeit.default_timer() - start_time:.2f}s)")

    # Write keep list if requested
    if keep_read_file:
        with open(keep_read_file, 'w') as f:
            for name in sorted(keep_set):
                f.write(f"{name}\n")

    # Phase 2: Write filtered BAM from original to_remap input
    print(f"Phase 2: Writing filtered BAM...")
    phase2_start = timeit.default_timer()
    kept_written = 0

    with pysam.AlignmentFile(to_remap_bam, "rb") as in_bam:
        header = in_bam.header
        with pysam.AlignmentFile(filt_out_bam, "wb", header=header, threads=threads) as out_bam:
            for rec in in_bam:
                if rec.query_name in keep_set:
                    out_bam.write(rec)
                    kept_written += 1

    print(f"Wrote {kept_written} reads to {filt_out_bam} (Phase 2: {timeit.default_timer() - phase2_start:.2f}s)")
    print(f"Total filtering time: {timeit.default_timer() - start_time:.2f}s")


def filt_remapped_reads(
    to_remap_bam: str,
    remapped_bam: str,
    filt_out_bam: str,
    keep_read_file: Optional[str] = None,
    use_rust: bool = True,
    threads: int = 1,
    same_locus_slop: int = 0,
) -> None:
    """Filter remapped reads using WASP algorithm.

    Args:
        to_remap_bam: Original BAM with reads to remap
        remapped_bam: Remapped BAM with swapped alleles
        filt_out_bam: Output filtered BAM
        keep_read_file: Optional file to write kept read names
        use_rust: Use Rust acceleration if available
        threads: Number of threads for BAM I/O
        same_locus_slop: Tolerance (bp) for same locus test (for indels)
    """
    rust_allowed = (
        use_rust
        and filter_bam_wasp is not None
        and os.environ.get("WASP2_DISABLE_RUST") != "1"
    )

    if not rust_allowed:
        print("⚠️  Rust extension not available, using Python fallback (slower)")
        _filter_bam_wasp_python(
            to_remap_bam,
            remapped_bam,
            filt_out_bam,
            keep_read_file=keep_read_file,
            threads=threads,
            same_locus_slop=same_locus_slop,
        )
        return

    # Rust path
    try:
        # Try with same_locus_slop (newer version)
        filter_bam_wasp(
            to_remap_bam,
            remapped_bam,
            filt_out_bam,
            keep_read_file=keep_read_file,
            threads=threads,
            same_locus_slop=same_locus_slop,
        )
    except TypeError:
        # Fall back to older version without same_locus_slop
        filter_bam_wasp(
            to_remap_bam,
            remapped_bam,
            filt_out_bam,
            keep_read_file=keep_read_file,
            threads=threads,
        )


def merge_filt_bam(
    keep_bam: str,
    remapped_filt_bam: str,
    out_bam: str,
    threads: int = 1
) -> None:
    """Merge filtered BAM files using samtools (faster than pysam).

    Both input BAMs are already coordinate-sorted, so samtools merge
    produces sorted output without needing an explicit sort step.

    Args:
        keep_bam: BAM with reads that didn't need remapping
        remapped_filt_bam: BAM with filtered remapped reads
        out_bam: Output merged BAM
        threads: Number of threads for samtools
    """
    start_time = timeit.default_timer()

    # Merge using samtools (faster than pysam, inputs are already sorted)
    subprocess.run(
        ["samtools", "merge", "-@", str(threads),
         "-f", "-o", out_bam, keep_bam, remapped_filt_bam],
        check=True)
    print(f"Merged BAM in {timeit.default_timer() - start_time:.2f} seconds")

    # Index the merged BAM (no sort needed - inputs were already sorted)
    start_index = timeit.default_timer()
    subprocess.run(
        ["samtools", "index", "-@", str(threads), out_bam],
        check=True)
    print(f"Indexed BAM in {timeit.default_timer() - start_index:.2f} seconds")
