#!/usr/bin/env python3
"""
Test Rust mapping implementation against Python baseline
"""
import sys
import os
import time
from pathlib import Path

# Add paths
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "rust" / "target" / "release"))
os.chdir(project_root)

try:
    import wasp2_rust
    print("✓ Rust module loaded successfully")
except ImportError as e:
    print(f"✗ Failed to load Rust module: {e}")
    sys.exit(1)

# Test files
bam_file = "baselines/mapping/to_remap.bam"
intersect_bed = "baselines/mapping/intersect.bed"
chrom = "chr10"

# Output files
rust_r1 = "/tmp/rust_mapping_r1.fq"
rust_r2 = "/tmp/rust_mapping_r2.fq"

print("\n" + "=" * 80)
print("RUST MAPPING TEST")
print("=" * 80)
print(f"BAM file: {bam_file}")
print(f"Intersect BED: {intersect_bed}")
print(f"Chromosome: {chrom}")
print()

# Test Rust implementation
print("Testing Rust implementation...")
start = time.time()
try:
    pairs_processed, haplotypes_generated = wasp2_rust.remap_chromosome(
        bam_file,
        intersect_bed,
        chrom,
        rust_r1,
        rust_r2,
        max_seqs=64
    )
    elapsed = time.time() - start

    print(f"✓ Success!")
    print(f"  Pairs processed: {pairs_processed}")
    print(f"  Haplotypes generated: {haplotypes_generated}")
    print(f"  Time: {elapsed:.4f}s")

    # Check output files
    if os.path.exists(rust_r1):
        r1_size = os.path.getsize(rust_r1)
        with open(rust_r1) as f:
            r1_lines = sum(1 for _ in f)
        print(f"  R1 output: {r1_size} bytes, {r1_lines} lines, {r1_lines//4} reads")
    else:
        print(f"  ✗ R1 output not found: {rust_r1}")

    if os.path.exists(rust_r2):
        r2_size = os.path.getsize(rust_r2)
        with open(rust_r2) as f:
            r2_lines = sum(1 for _ in f)
        print(f"  R2 output: {r2_size} bytes, {r2_lines} lines, {r2_lines//4} reads")
    else:
        print(f"  ✗ R2 output not found: {rust_r2}")

    # Show sample output
    if os.path.exists(rust_r1):
        print(f"\n  Sample R1 output (first 8 lines):")
        with open(rust_r1) as f:
            for i, line in enumerate(f):
                if i >= 8:
                    break
                print(f"    {line.rstrip()}")

except Exception as e:
    elapsed = time.time() - start
    print(f"✗ Failed: {e}")
    print(f"  Time: {elapsed:.4f}s")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 80)
print("TEST COMPLETE")
print("=" * 80)
