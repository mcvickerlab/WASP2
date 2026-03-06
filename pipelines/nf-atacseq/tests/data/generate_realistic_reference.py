#!/usr/bin/env python3
"""
Generate a realistic ~20kb non-repetitive reference sequence for ATAC-seq testing.

Properties:
  - ~42% GC content (human-like)
  - No homopolymer runs > 5bp
  - High k-mer uniqueness (>99% unique 20-mers)
  - Deterministic (seeded RNG)

Output: chr_test.fa with contig name 'chr_test'
"""

import random
import sys

SEED = 42
LENGTH = 20000
CONTIG = "chr_test"
LINE_WIDTH = 60
MAX_HOMOPOLYMER = 5

# Target base frequencies for ~42% GC
# A=29%, T=29%, G=21%, C=21%
BASES = "ATGC"
WEIGHTS = [0.29, 0.29, 0.21, 0.21]


def generate_sequence(length, seed=SEED):
    """Generate a non-repetitive sequence with controlled GC content."""
    rng = random.Random(seed)

    seq = []
    homopolymer_count = 0
    last_base = None

    for _ in range(length):
        # Pick a base using weighted random
        base = rng.choices(BASES, weights=WEIGHTS, k=1)[0]

        # Prevent long homopolymers
        if base == last_base:
            homopolymer_count += 1
            if homopolymer_count >= MAX_HOMOPOLYMER:
                # Force a different base
                alternatives = [b for b in BASES if b != base]
                alt_weights = [WEIGHTS[BASES.index(b)] for b in alternatives]
                total = sum(alt_weights)
                alt_weights = [w / total for w in alt_weights]
                base = rng.choices(alternatives, weights=alt_weights, k=1)[0]
                homopolymer_count = 1
        else:
            homopolymer_count = 1

        seq.append(base)
        last_base = base

    return "".join(seq)


def validate_sequence(seq):
    """Validate sequence properties."""
    gc = sum(1 for b in seq if b in "GC") / len(seq)

    # Check k-mer uniqueness
    kmers_20 = set()
    for i in range(len(seq) - 19):
        kmers_20.add(seq[i : i + 20])
    unique_20 = len(kmers_20)
    total_20 = len(seq) - 19
    uniqueness = unique_20 / total_20

    # Check max homopolymer
    max_hp = 1
    current_hp = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_hp += 1
            max_hp = max(max_hp, current_hp)
        else:
            current_hp = 1

    return {
        "length": len(seq),
        "gc_content": gc,
        "unique_20mers": unique_20,
        "total_20mers": total_20,
        "uniqueness_pct": uniqueness * 100,
        "max_homopolymer": max_hp,
    }


def write_fasta(seq, contig, filepath, line_width=LINE_WIDTH):
    """Write sequence as FASTA."""
    with open(filepath, "w") as f:
        f.write(f">{contig}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i : i + line_width] + "\n")


def main():
    output = sys.argv[1] if len(sys.argv) > 1 else "chr_test.fa"

    print(f"Generating {LENGTH}bp non-repetitive reference sequence...")
    seq = generate_sequence(LENGTH)

    stats = validate_sequence(seq)
    print(f"  Length: {stats['length']}bp")
    print(f"  GC content: {stats['gc_content']:.1%}")
    print(f"  Unique 20-mers: {stats['unique_20mers']}/{stats['total_20mers']} ({stats['uniqueness_pct']:.1f}%)")
    print(f"  Max homopolymer: {stats['max_homopolymer']}bp")

    # Validate
    assert stats["gc_content"] > 0.38 and stats["gc_content"] < 0.46, f"GC content out of range: {stats['gc_content']}"
    assert stats["uniqueness_pct"] > 99.0, f"Uniqueness too low: {stats['uniqueness_pct']}"
    assert stats["max_homopolymer"] <= MAX_HOMOPOLYMER, f"Homopolymer too long: {stats['max_homopolymer']}"

    write_fasta(seq, CONTIG, output)
    print(f"  Wrote {output}")


if __name__ == "__main__":
    main()
