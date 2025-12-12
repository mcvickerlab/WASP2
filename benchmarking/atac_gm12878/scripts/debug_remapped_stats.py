#!/usr/bin/env python3
"""Compare remapped BAM statistics between Python and Rust."""

import pysam
from collections import Counter, defaultdict
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


def analyze_remapped_bam(bam_path):
    """Analyze remapped BAM and return statistics."""
    stats = {
        'total_records': 0,
        'proper_pairs': 0,
        'secondary': 0,
        'supplementary': 0,
        'unmapped': 0,
        'unique_wasp_names': set(),
        'unique_orig_names': set(),
        'total_values': Counter(),  # histogram of 'total' values from WASP names
        'pairs_per_orig': defaultdict(int),  # how many pairs per original read name
    }

    read_dict = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            stats['total_records'] += 1

            if read.is_unmapped:
                stats['unmapped'] += 1
            if read.is_proper_pair:
                stats['proper_pairs'] += 1
            if read.is_secondary:
                stats['secondary'] += 1
            if read.is_supplementary:
                stats['supplementary'] += 1

            # Filter for counting (match filter logic)
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name
            parsed = parse_wasp_name(qname)
            if parsed is None:
                continue

            orig_name, pos1, pos2, seq_idx, total = parsed
            stats['unique_wasp_names'].add(qname)
            stats['unique_orig_names'].add(orig_name)

            # Count pairs
            if qname not in read_dict:
                read_dict[qname] = read
            else:
                read_dict.pop(qname)
                stats['total_values'][total] += 1
                stats['pairs_per_orig'][orig_name] += 1

    stats['incomplete_pairs'] = len(read_dict)
    stats['unique_wasp_names'] = len(stats['unique_wasp_names'])
    stats['unique_orig_names'] = len(stats['unique_orig_names'])

    return stats


def main():
    python_bam = "results/wasp2python_snp_fixed_2025-12-07_13-12-00/remapped.bam"
    rust_bam = "results/wasp2rust_snp_fixed_2025-12-09_22-40-17/remapped.bam"

    print("=" * 80)
    print("Remapped BAM Comparison")
    print("=" * 80)

    print(f"\nAnalyzing Python remapped BAM: {python_bam}")
    py_stats = analyze_remapped_bam(python_bam)

    print(f"\nAnalyzing Rust remapped BAM: {rust_bam}")
    rust_stats = analyze_remapped_bam(rust_bam)

    print("\n" + "=" * 80)
    print("COMPARISON")
    print("=" * 80)

    metrics = [
        ('total_records', 'Total records'),
        ('proper_pairs', 'Proper pairs'),
        ('secondary', 'Secondary'),
        ('supplementary', 'Supplementary'),
        ('unmapped', 'Unmapped'),
        ('unique_wasp_names', 'Unique WASP names'),
        ('unique_orig_names', 'Unique original names'),
        ('incomplete_pairs', 'Incomplete pairs'),
    ]

    print(f"\n{'Metric':<25} {'Python':>15} {'Rust':>15} {'Diff':>15}")
    print("-" * 70)
    for key, label in metrics:
        py_val = py_stats[key]
        rust_val = rust_stats[key]
        diff = rust_val - py_val
        print(f"{label:<25} {py_val:>15,} {rust_val:>15,} {diff:>+15,}")

    print("\n--- Total values histogram (pairs per haplotype count) ---")
    print(f"{'Total':<10} {'Python pairs':>15} {'Rust pairs':>15} {'Diff':>15}")
    print("-" * 55)

    all_totals = set(py_stats['total_values'].keys()) | set(rust_stats['total_values'].keys())
    for t in sorted(all_totals)[:10]:  # Show first 10
        py_val = py_stats['total_values'].get(t, 0)
        rust_val = rust_stats['total_values'].get(t, 0)
        diff = rust_val - py_val
        print(f"{t:<10} {py_val:>15,} {rust_val:>15,} {diff:>+15,}")

    # Check if reads with different 'total' counts are handled differently
    print("\n--- Pairs per original read name histogram ---")
    py_pairs_hist = Counter(py_stats['pairs_per_orig'].values())
    rust_pairs_hist = Counter(rust_stats['pairs_per_orig'].values())

    print(f"{'Pairs':<10} {'Python orig':>15} {'Rust orig':>15} {'Diff':>15}")
    print("-" * 55)
    all_pairs = set(py_pairs_hist.keys()) | set(rust_pairs_hist.keys())
    for p in sorted(all_pairs)[:10]:
        py_val = py_pairs_hist.get(p, 0)
        rust_val = rust_pairs_hist.get(p, 0)
        diff = rust_val - py_val
        print(f"{p:<10} {py_val:>15,} {rust_val:>15,} {diff:>+15,}")


if __name__ == "__main__":
    main()
