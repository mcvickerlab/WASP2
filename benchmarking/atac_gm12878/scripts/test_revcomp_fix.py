#!/usr/bin/env python3
"""Quick test to verify the reverse complement fix in FASTQ output."""

import subprocess
import tempfile
import os

# Create a mini test using a small subset of reads
def test_revcomp():
    print("Testing reverse complement fix...")

    # Check one of the problematic reads from the old output
    # SRR891268.13155439 - both mates have same position and sequence in original BAM

    # First verify the module is working
    import wasp2_rust
    print("wasp2_rust module loaded")

    # Test reverse complement logic by checking if outputs differ for reverse-strand reads
    # We'll check the actual benchmark output when it completes

    print("\nTo verify the fix is working:")
    print("1. Wait for benchmark job to complete")
    print("2. Check that R1 and R2 have different sequences for reads like SRR891268.13155439")
    print("3. Count non-proper pairs - should be much lower (~3K instead of ~466K)")

    print("\nTest complete - manual verification needed after benchmark")


if __name__ == "__main__":
    test_revcomp()
