#!/usr/bin/env python3
"""
Extract SNV read counts from benchmark JSON files for Figure 1C.
This script reads the benchmark results and generates a TSV file for plotting.
"""

import json
from pathlib import Path
import sys

# Define paths
RESULTS_DIR = Path(__file__).parent.parent / "results"
OUTPUT_FILE = Path(__file__).parent.parent.parent.parent / "paper/figure1/data/snv_indel_read_counts.tsv"

# Define which benchmark directories to use
BENCHMARKS = {
    "wasp1": "wasp1_snp_FIXED_2025-12-07_12-27-45",
    "wasp2python": "wasp2python_snp_DEV_MT_2025-12-07_15-52-30",
    "wasp2rust_snp": "wasp2rust_snp_fixed_2025-12-06_22-38-30",
    "wasp2rust_indel": "wasp2rust_indel_fixed_2025-12-07_12-51-31",
}

def extract_snv_counts():
    """Extract SNV read counts from all benchmark JSON files."""
    results = []

    for pipeline, dirname in BENCHMARKS.items():
        benchmark_dir = RESULTS_DIR / dirname
        json_file = benchmark_dir / "benchmark_results.json"

        if not json_file.exists():
            print(f"WARNING: {json_file} not found, skipping {pipeline}")
            continue

        with open(json_file) as f:
            data = json.load(f)

        # Extract SNV counts - field names vary slightly between pipelines
        snv_pre = data.get("snv_reads_pre")
        snv_post = data.get("snv_reads_post")

        # Calculate pass rate if not present
        if snv_pre and snv_post:
            pass_rate = (snv_post / snv_pre) * 100
        else:
            pass_rate = data.get("snv_read_pass_rate") or data.get("snv_pass_rate_percent")

        if snv_pre is None or snv_post is None:
            print(f"WARNING: Missing SNV counts in {json_file}")
            print(f"  Available keys: {list(data.keys())}")
            continue

        results.append({
            "pipeline": pipeline,
            "variant_type": "SNV",
            "reads_pre": snv_pre,
            "reads_post": snv_post,
            "pass_rate": round(pass_rate, 2)
        })

        print(f"Extracted from {pipeline}:")
        print(f"  JSON file: {json_file}")
        print(f"  SNV pre:  {snv_pre:,}")
        print(f"  SNV post: {snv_post:,}")
        print(f"  Pass rate: {pass_rate:.2f}%")
        print()

    return results

def write_tsv(results, output_file):
    """Write results to TSV file."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        # Header
        f.write("pipeline\tvariant_type\treads_pre\treads_post\tpass_rate\n")

        # Data rows
        for row in results:
            f.write(f"{row['pipeline']}\t{row['variant_type']}\t{row['reads_pre']}\t{row['reads_post']}\t{row['pass_rate']}\n")

    print(f"Written to: {output_file}")

def main():
    print("=" * 60)
    print("Extracting SNV Read Counts for Figure 1C")
    print("=" * 60)
    print()

    results = extract_snv_counts()

    if not results:
        print("ERROR: No data extracted!")
        sys.exit(1)

    print("=" * 60)
    print("Summary")
    print("=" * 60)
    for r in results:
        print(f"{r['pipeline']:20s} | {r['reads_pre']:>10,} -> {r['reads_post']:>10,} | {r['pass_rate']:>6.2f}%")

    print()
    write_tsv(results, OUTPUT_FILE)

if __name__ == "__main__":
    main()
