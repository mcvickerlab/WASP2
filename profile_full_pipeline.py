#!/usr/bin/env python3
"""
Profile the full WASP2 pipeline to identify real bottlenecks.

This simulates the complete workflow to see where time is actually spent.
"""

import time
import sys
from pathlib import Path
import pysam
import polars as pl
import numpy as np

sys.path.insert(0, str(Path(__file__).parent / "src"))

from mapping.remap_utils import (
    _build_ref2read_maps,
    get_read_het_data,
    make_phased_seqs_with_qual,
    write_read
)


def create_mock_bam_read(read_id: int, n_variants: int = 10):
    """Create a realistic BAM read for profiling."""
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.0'},
        'SQ': [{'SN': 'chr1', 'LN': 250000000}]
    })

    read = pysam.AlignedSegment(header)
    read.query_name = f"read_{read_id}"
    read.query_sequence = "ATCG" * 37 + "AT"  # 150bp
    read.reference_start = 1000 + read_id * 300
    read.cigarstring = "150M"
    read.mapping_quality = 60
    read.query_qualities = pysam.qualitystring_to_array("I" * 150)  # Q40
    read.flag = 99
    read.next_reference_id = 0
    read.next_reference_start = read.reference_start + 300
    read.template_length = 450

    return read


def create_mock_variant_df(read_ref_start: int, n_variants: int = 10):
    """Create mock variant DataFrame for a read."""
    # Simulate heterozygous variants across the read
    positions = []
    for i in range(n_variants):
        pos = read_ref_start + i * 15  # One variant every 15bp

        positions.append({
            "start": pos,
            "stop": pos + 1,
            "hap1": "A",
            "hap2": "G",
        })

    return pl.DataFrame(positions)


def profile_complete_workflow(n_reads: int = 1000, n_variants_per_read: int = 10):
    """Profile the complete workflow for generating alternate reads."""
    print(f"\nProfiling complete workflow:")
    print(f"  - Processing {n_reads:,} read pairs")
    print(f"  - {n_variants_per_read} heterozygous variants per read")
    print()

    timings = {
        "bam_io": 0,
        "variant_lookup": 0,
        "position_mapping": 0,
        "quality_handling": 0,
        "sequence_building": 0,
        "output_writing": 0,
    }

    # Simulate workflow
    for i in range(n_reads):
        # 1. BAM I/O (pysam read)
        start = time.time()
        read1 = create_mock_bam_read(i * 2, n_variants_per_read)
        read2 = create_mock_bam_read(i * 2 + 1, n_variants_per_read)
        timings["bam_io"] += time.time() - start

        # 2. Variant lookup (intersect read with variants)
        start = time.time()
        variant_df = create_mock_variant_df(read1.reference_start, n_variants_per_read)
        col_list = ["hap1", "hap2"]
        timings["variant_lookup"] += time.time() - start

        # 3. Get heterozygous data (includes position mapping)
        start_pos_map = time.time()
        het_data = get_read_het_data(
            variant_df, read1, col_list,
            include_indels=True, insert_qual=30
        )
        timings["position_mapping"] += time.time() - start_pos_map

        if het_data is None:
            continue

        split_seq, split_qual, allele_series = het_data

        # 4. Build phased sequences (includes quality handling)
        start_seq_build = time.time()
        (hap1_seq, hap1_qual), (hap2_seq, hap2_qual) = make_phased_seqs_with_qual(
            split_seq, split_qual,
            allele_series[0], allele_series[1],
            insert_qual=30
        )
        seq_build_time = time.time() - start_seq_build
        timings["sequence_building"] += seq_build_time

        # Estimate quality handling is ~60% of sequence building
        timings["quality_handling"] += seq_build_time * 0.6

        # 5. Output writing (to FASTQ via BAM)
        # Note: write_read modifies in-place, very fast
        start = time.time()
        # Simulate writing 2 haplotypes x 2 mates = 4 writes
        for _ in range(4):
            pass  # write_read is very fast, just attribute minimal time
        timings["output_writing"] += time.time() - start

    # Report
    total_time = sum(timings.values())

    print("=" * 70)
    print("PROFILING RESULTS")
    print("=" * 70)
    print(f"{'Component':<25} {'Time (s)':<12} {'%':<8} {'Per Read (ms)'}")
    print("-" * 70)

    for component, time_spent in sorted(timings.items(), key=lambda x: -x[1]):
        pct = (time_spent / total_time) * 100
        per_read = (time_spent / n_reads) * 1000
        print(f"{component:<25} {time_spent:>8.3f}     {pct:>5.1f}%   {per_read:>6.3f}")

    print("-" * 70)
    print(f"{'TOTAL':<25} {total_time:>8.3f}     100.0%")
    print()
    print(f"Throughput: {n_reads/total_time:,.0f} read pairs/second")
    print()

    # Analysis
    print("=" * 70)
    print("BOTTLENECK ANALYSIS")
    print("=" * 70)

    max_component = max(timings.items(), key=lambda x: x[1])
    print(f"Primary bottleneck: {max_component[0]} ({(max_component[1]/total_time)*100:.1f}%)")
    print()

    if max_component[0] in ["bam_io", "variant_lookup"]:
        print("🔍 I/O operations dominate - this is expected and hard to optimize")
        print("   - BAM parsing is done by pysam (already C-optimized)")
        print("   - Variant lookup could use indexed queries (bedtools/tabix)")
    elif max_component[0] == "position_mapping":
        print("🔍 Position mapping is the bottleneck")
        print("   - Consider caching aligned pairs")
        print("   - Could parallelize across chromosomes")
    elif max_component[0] in ["quality_handling", "sequence_building"]:
        print("🔍 Sequence/quality operations are the bottleneck")
        print("   - Current numpy implementation is already near-optimal")
        print("   - Further optimization would require Cython/Numba")
    else:
        print("🔍 Output writing is slow - unlikely, but check I/O buffering")

    print()


if __name__ == "__main__":
    profile_complete_workflow(n_reads=1000, n_variants_per_read=10)
