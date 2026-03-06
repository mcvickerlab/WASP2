#!/usr/bin/env python3
"""Generate a realistic ~20kb synthetic reference genome for WASP2 tests.

Properties:
- Single contig named 'chr_test'
- ~40-45% GC content (realistic for mammalian genomes)
- No homopolymer runs longer than 5bp
- High sequence complexity for unique k-mer mapping
- Deterministic output (fixed seed)
"""

import random
import sys


def generate_reference(length=19800, seed=12345, max_homopolymer=5, line_width=60):
    """Generate a random reference sequence with realistic properties."""
    rng = random.Random(seed)

    # Weighted nucleotide frequencies for ~42% GC content
    # A=29%, T=29%, G=21%, C=21%
    bases = ['A', 'T', 'G', 'C']
    weights = [0.29, 0.29, 0.21, 0.21]

    sequence = []
    run_count = 0
    last_base = None

    for _ in range(length):
        # Pick a base with the desired frequency distribution
        base = rng.choices(bases, weights=weights, k=1)[0]

        # Enforce max homopolymer constraint
        if base == last_base:
            run_count += 1
            if run_count >= max_homopolymer:
                # Force a different base
                other_bases = [b for b in bases if b != base]
                other_weights = [w for b, w in zip(bases, weights) if b != base]
                total = sum(other_weights)
                other_weights = [w / total for w in other_weights]
                base = rng.choices(other_bases, weights=other_weights, k=1)[0]
                run_count = 1
        else:
            run_count = 1

        last_base = base
        sequence.append(base)

    seq_str = ''.join(sequence)

    # Verify properties
    gc_count = seq_str.count('G') + seq_str.count('C')
    gc_pct = gc_count / len(seq_str) * 100

    # Check max homopolymer
    max_run = 0
    current_run = 1
    for i in range(1, len(seq_str)):
        if seq_str[i] == seq_str[i-1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1

    print(f"Reference stats:", file=sys.stderr)
    print(f"  Length: {len(seq_str)} bp", file=sys.stderr)
    print(f"  GC content: {gc_pct:.1f}%", file=sys.stderr)
    print(f"  Max homopolymer: {max_run} bp", file=sys.stderr)

    # Write FASTA
    print(">chr_test")
    for i in range(0, len(seq_str), line_width):
        print(seq_str[i:i+line_width])

    return seq_str


def extract_bases_at_positions(seq_str, positions):
    """Print the base at each 1-based position (for VCF REF allele verification)."""
    print("\nBases at SNP positions (1-based):", file=sys.stderr)
    for pos in sorted(positions):
        if 1 <= pos <= len(seq_str):
            base = seq_str[pos - 1]  # Convert to 0-based
            print(f"  pos {pos}: {base}", file=sys.stderr)


if __name__ == '__main__':
    # VCF SNP positions from the test data
    snp_positions = [
        750, 1200, 2800, 3200, 5000,
        6000, 6100, 6700, 6800, 7400, 7500,
        8100, 8200, 8800, 8900,
        10800, 11200, 12800, 13200, 15000,
        16000, 16100, 16700, 16800, 17400, 17500,
        18100, 18200, 18800, 18900,
    ]

    seq = generate_reference()
    extract_bases_at_positions(seq, snp_positions)
