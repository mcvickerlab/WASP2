#!/usr/bin/env python3
"""Debug script to trace exact output from filters."""

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


def python_filter_logic(remapped_bam_path):
    """Replicate Python filt_remapped_reads logic exactly, return keep_set."""
    pos_dict = {}
    total_dict = {}
    keep_set = set()

    read_dict = {}

    with pysam.AlignmentFile(remapped_bam_path, "rb") as bam:
        for read in bam.fetch():
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name

            if qname not in read_dict:
                read_dict[qname] = read
                continue

            if read.is_read1:
                read1 = read
                read2 = read_dict.pop(qname)
            else:
                read1 = read_dict.pop(qname)
                read2 = read

            parsed = parse_wasp_name(qname)
            if parsed is None:
                continue

            orig_name, pos1, pos2, seq_idx, total = parsed

            if orig_name not in pos_dict:
                pos_dict[orig_name] = (pos1, pos2)
                total_dict[orig_name] = total
                keep_set.add(orig_name)
            elif orig_name not in keep_set:
                continue

            total_dict[orig_name] -= 1

            actual_pos = (read1.reference_start, read1.next_reference_start)
            expected_pos = pos_dict[orig_name]

            if actual_pos != expected_pos:
                keep_set.remove(orig_name)
                total_dict.pop(orig_name)
            elif total_dict[orig_name] == 0:
                total_dict.pop(orig_name)
                pos_dict.pop(orig_name)

    missing_count_set = set(total_dict.keys())
    keep_set = keep_set - missing_count_set

    return keep_set


def count_reads_in_to_remap(to_remap_bam, keep_set):
    """Count how many reads from to_remap BAM match the keep_set."""
    matched = 0
    total = 0

    with pysam.AlignmentFile(to_remap_bam, "rb") as bam:
        for read in bam.fetch():
            total += 1
            if read.query_name in keep_set:
                matched += 1

    return matched, total


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--remapped", required=True)
    parser.add_argument("--to_remap", required=True)
    args = parser.parse_args()

    print(f"Remapped BAM: {args.remapped}")
    print(f"To-remap BAM: {args.to_remap}")
    print()

    print("Building keep_set from remapped BAM...")
    keep_set = python_filter_logic(args.remapped)
    print(f"Keep set has {len(keep_set)} original read names")
    print()

    print("Counting matching reads in to_remap BAM...")
    matched, total = count_reads_in_to_remap(args.to_remap, keep_set)
    print(f"To-remap BAM total reads: {total}")
    print(f"Reads matching keep_set: {matched}")
    print(f"Match rate: {matched/total*100:.2f}%")


if __name__ == "__main__":
    main()
