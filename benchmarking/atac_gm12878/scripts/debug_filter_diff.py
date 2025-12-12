#!/usr/bin/env python3
"""Debug script to identify why Rust filter is stricter than Python filter.

Compare position checking logic between:
- Python: (read1.reference_start, read1.next_reference_start) == (pos1, pos2)
- Rust: checks both orderings, uses second arriving record
"""

import pysam
from collections import defaultdict
import sys

def parse_wasp_name(qname):
    """Parse WASP-encoded read name."""
    if "_WASP_" not in qname:
        return None
    parts = qname.split("_WASP_")
    orig_name = parts[0]
    data = parts[1].split("_")
    pos1 = int(data[0])
    pos2 = int(data[1])
    seq_idx = int(data[2])
    total = int(data[3])
    return orig_name, pos1, pos2, seq_idx, total


def python_filter_logic(remapped_bam_path, verbose=False):
    """Replicate Python filt_remapped_reads logic exactly."""
    pos_dict = {}
    total_dict = {}
    keep_set = set()
    removed_details = []

    # Track read pairs (mimic paired_read_gen)
    read_dict = {}

    with pysam.AlignmentFile(remapped_bam_path, "rb") as bam:
        for read in bam.fetch():
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name

            if qname not in read_dict:
                read_dict[qname] = read
                continue

            # Pair complete - determine read1 and read2
            if read.is_read1:
                read1 = read
                read2 = read_dict.pop(qname)
            else:
                read1 = read_dict.pop(qname)
                read2 = read

            # Parse WASP name
            parsed = parse_wasp_name(qname)
            if parsed is None:
                continue

            orig_name, pos1, pos2, seq_idx, total = parsed

            # Initialize tracking
            if orig_name not in pos_dict:
                pos_dict[orig_name] = (pos1, pos2)
                total_dict[orig_name] = total
                keep_set.add(orig_name)
            elif orig_name not in keep_set:
                continue

            # Count down
            total_dict[orig_name] -= 1

            # Python comparison: (read1.reference_start, read1.next_reference_start)
            actual_pos = (read1.reference_start, read1.next_reference_start)
            expected_pos = pos_dict[orig_name]

            if actual_pos != expected_pos:
                keep_set.remove(orig_name)
                total_dict.pop(orig_name)
                if verbose and len(removed_details) < 20:
                    removed_details.append({
                        'name': orig_name,
                        'wasp_qname': qname,
                        'expected': expected_pos,
                        'actual': actual_pos,
                        'read1_pos': read1.reference_start,
                        'read1_mpos': read1.next_reference_start,
                        'read2_pos': read2.reference_start,
                        'read2_mpos': read2.next_reference_start,
                        'reason': 'position_mismatch'
                    })
            elif total_dict[orig_name] == 0:
                total_dict.pop(orig_name)
                pos_dict.pop(orig_name)

    # Remove reads with missing counts
    missing_count_set = set(total_dict.keys())
    for name in missing_count_set:
        if verbose and len(removed_details) < 50:
            removed_details.append({
                'name': name,
                'remaining': total_dict[name],
                'reason': 'missing_count'
            })
    keep_set = keep_set - missing_count_set

    return keep_set, removed_details


def rust_filter_logic(remapped_bam_path, verbose=False):
    """Replicate Rust filter_bam_wasp logic exactly."""
    pos_map = {}
    remaining = {}
    keep_set = set()
    removed_details = []

    # Buffer for incomplete pairs
    read_buffer = {}

    with pysam.AlignmentFile(remapped_bam_path, "rb") as bam:
        for rec in bam.fetch():
            if rec.is_unmapped or not rec.is_proper_pair or rec.is_secondary or rec.is_supplementary:
                continue

            qname = rec.query_name

            parsed = parse_wasp_name(qname)
            if parsed is None:
                continue

            orig_name, pos1, pos2, seq_idx, total = parsed

            # Buffer until both mates arrive (using full qname including WASP suffix)
            if qname not in read_buffer:
                read_buffer[qname] = {
                    'pos': rec.reference_start,
                    'mpos': rec.next_reference_start,
                    'is_read1': rec.is_read1
                }
                continue

            # Second mate arrived
            first_read = read_buffer.pop(qname)

            # Initialize tracking
            if orig_name not in pos_map:
                pos_map[orig_name] = (pos1, pos2)
                remaining[orig_name] = total
                keep_set.add(orig_name)
            elif orig_name not in keep_set:
                continue

            # Count down
            remaining[orig_name] -= 1

            # Rust comparison: uses second arriving record's positions
            rec_pos = rec.reference_start
            mate_pos = rec.next_reference_start
            expect_pos, expect_mate = pos_map[orig_name]

            # Check both orderings
            matches = (
                (rec_pos == expect_pos and mate_pos == expect_mate) or
                (rec_pos == expect_mate and mate_pos == expect_pos)
            )

            if not matches:
                keep_set.remove(orig_name)
                if verbose and len(removed_details) < 20:
                    removed_details.append({
                        'name': orig_name,
                        'wasp_qname': qname,
                        'expected': (expect_pos, expect_mate),
                        'actual_rec': (rec_pos, mate_pos),
                        'first_read': first_read,
                        'second_is_read1': rec.is_read1,
                        'reason': 'position_mismatch'
                    })
                continue

            # Drop bookkeeping if all pairs seen
            if remaining[orig_name] <= 0:
                remaining.pop(orig_name)
                pos_map.pop(orig_name)

    # Remove reads with missing counts
    for name in list(remaining.keys()):
        keep_set.discard(name)
        if verbose and len(removed_details) < 50:
            removed_details.append({
                'name': name,
                'remaining': remaining[name],
                'reason': 'missing_count'
            })

    return keep_set, removed_details


def main():
    import sys

    # Compare both benchmarks
    rust_remapped = "results/wasp2rust_snp_fixed_2025-12-09_22-40-17/remapped.bam"
    python_remapped = "results/wasp2python_snp_fixed_2025-12-07_13-12-00/remapped.bam"

    # Use command line arg to select which one
    if len(sys.argv) > 1:
        if sys.argv[1] == "python":
            remapped_bam = python_remapped
            print("Using PYTHON benchmark remapped.bam")
        elif sys.argv[1] == "rust":
            remapped_bam = rust_remapped
            print("Using RUST benchmark remapped.bam")
        else:
            remapped_bam = sys.argv[1]
    else:
        remapped_bam = rust_remapped
        print("Using RUST benchmark remapped.bam (default)")

    print("=" * 80)
    print("Filter Logic Comparison: Python vs Rust")
    print("=" * 80)

    print(f"\nReading remapped BAM: {remapped_bam}")

    print("\n--- Running Python filter logic ---")
    python_keep, python_removed = python_filter_logic(remapped_bam, verbose=True)
    print(f"Python keeps: {len(python_keep)} read pairs")

    print("\n--- Running Rust filter logic ---")
    rust_keep, rust_removed = rust_filter_logic(remapped_bam, verbose=True)
    print(f"Rust keeps: {len(rust_keep)} read pairs")

    print("\n" + "=" * 80)
    print("COMPARISON")
    print("=" * 80)

    only_python = python_keep - rust_keep
    only_rust = rust_keep - python_keep

    print(f"Both keep: {len(python_keep & rust_keep)}")
    print(f"Only Python keeps: {len(only_python)}")
    print(f"Only Rust keeps: {len(only_rust)}")
    print(f"Difference (Python - Rust): {len(python_keep) - len(rust_keep)}")

    # Show some examples of reads only Python keeps
    if only_python:
        print("\n--- Examples: Reads Python keeps but Rust rejects ---")
        examples = list(only_python)[:10]
        for name in examples:
            print(f"  {name}")

    # Show removed details
    if python_removed:
        print("\n--- Python removal details (first 10) ---")
        for d in python_removed[:10]:
            if d['reason'] == 'position_mismatch':
                print(f"  {d['name']}: expected={d['expected']}, actual={d['actual']}")
                print(f"    read1: pos={d['read1_pos']}, mpos={d['read1_mpos']}")
                print(f"    read2: pos={d['read2_pos']}, mpos={d['read2_mpos']}")

    if rust_removed:
        print("\n--- Rust removal details (first 10) ---")
        for d in rust_removed[:10]:
            if d['reason'] == 'position_mismatch':
                print(f"  {d['name']}: expected={d['expected']}, actual_rec={d['actual_rec']}")


if __name__ == "__main__":
    main()
