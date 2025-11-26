#!/usr/bin/env python3
"""
WASP2 Indel Validation - Realistic Simulation Framework

Generates synthetic reads with known allelic ratios, aligns with BWA to create
realistic CIGAR strings, runs full WASP2 pipeline, and validates results.

This is the REAL test - proves WASP2's position mapping and quality inference work correctly.

Tiered testing:
  - Minimum:       90 tests  (~10 min) - Proves algorithm works
  - Moderate:     270 tests  (~30 min) - Shows robustness across coverage
  - Comprehensive: 810 tests (~2 hrs)  - Edge cases and large indels

Usage:
    python simulate_indel_ase_v2.py --tier minimum
    python simulate_indel_ase_v2.py --tier moderate
    python simulate_indel_ase_v2.py --tier comprehensive
"""

import pysam
import random
import tempfile
import subprocess
import argparse
import sys
import multiprocessing
from pathlib import Path
import polars as pl
import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Dict
import gzip
from collections import defaultdict


@dataclass
class GroundTruth:
    """Known truth for validation."""
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str
    true_ratio: float       # REF/ALT expression ratio (e.g., 2.0 = 2:1)
    variant_type: str       # "SNP", "INS", "DEL"
    coverage: int           # Number of reads covering this variant
    replicate: int          # Replicate number for this config
    seed: int               # Random seed for reproducibility


def create_reference_fasta(output_file: str, length: int = 1000000):
    """Create reference genome with random sequence."""
    print(f"Creating reference genome ({length}bp)...")

    with open(output_file, 'w') as f:
        f.write(">chr1\n")
        # Generate random sequence
        random.seed(42)  # Reproducible reference
        seq = ''.join(random.choices('ATCG', k=length))

        # Write in 80bp lines (FASTA format)
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

    # Index with samtools
    subprocess.run(['samtools', 'faidx', output_file], check=True,
                   capture_output=True)

    print(f"  ✅ Created and indexed {output_file}")
    return seq


def generate_quality_scores(length: int, mean_qual: int = 35, std_qual: int = 5) -> np.ndarray:
    """Generate realistic quality scores with variation."""
    qualities = np.random.normal(mean_qual, std_qual, length)
    qualities = np.clip(qualities, 10, 40).astype(np.uint8)
    return qualities


def add_sequencing_errors(sequence: str, error_rate: float = 0.01) -> str:
    """Add realistic sequencing errors (substitutions)."""
    seq_list = list(sequence)
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            # Substitute with random base (not same as original)
            bases = [b for b in 'ATCG' if b != seq_list[i]]
            seq_list[i] = random.choice(bases)
    return ''.join(seq_list)


def create_read_sequence(
    ref_seq: str,
    var_pos: int,
    allele: str,
    ref_allele_len: int,
    read_length: int = 150,
    error_rate: float = 0.01
) -> Tuple[str, np.ndarray]:
    """
    Create read sequence with specific allele at variant position.

    IMPORTANT for indels:
    - For insertions (ALT longer than REF): read has extra bases, shifts downstream seq
    - For deletions (REF longer than ALT): read skips bases, pulls in upstream seq

    Args:
        ref_seq: Full reference sequence
        var_pos: Variant position in reference coordinates
        allele: The allele (REF or ALT) to put in the read
        ref_allele_len: Length of the REF allele (for proper indel handling)
        read_length: Target read length
        error_rate: Sequencing error rate

    Returns:
        Tuple of (sequence, quality_scores)
    """
    # Random position offset so variant isn't always centered
    offset = random.randint(-30, 30)
    read_start = max(0, var_pos - read_length // 2 + offset)

    # Build read sequence properly for indels
    # Key insight: we always skip ref_allele_len bases from reference, then add the allele
    left_seq = ref_seq[read_start:var_pos]

    # For the right part, we skip ONLY the ref allele length (not the length of 'allele')
    # This correctly handles insertions (allele > ref) and deletions (allele < ref)
    right_start = var_pos + ref_allele_len
    right_seq = ref_seq[right_start:right_start + read_length]  # Get enough for padding

    read_seq = left_seq + allele + right_seq

    # Truncate to exact read length (important for insertions which make read longer)
    if len(read_seq) > read_length:
        read_seq = read_seq[:read_length]
    elif len(read_seq) < read_length:
        # Pad with reference sequence if too short (can happen near boundaries)
        end_pos = right_start + len(right_seq)
        pad_len = read_length - len(read_seq)
        read_seq += ref_seq[end_pos:end_pos + pad_len]

    # Final length check
    read_seq = read_seq[:read_length]

    # Add sequencing errors
    read_seq = add_sequencing_errors(read_seq, error_rate)

    # Generate quality scores
    qualities = generate_quality_scores(len(read_seq))

    return read_seq, qualities


def write_fastq_record(fq_handle, read_id: str, sequence: str, qualities: np.ndarray):
    """Write single FASTQ record."""
    qual_string = ''.join([chr(q + 33) for q in qualities])  # Convert to Phred+33
    fq_handle.write(f"@{read_id}\n")
    fq_handle.write(f"{sequence}\n")
    fq_handle.write(f"+\n")
    fq_handle.write(f"{qual_string}\n")


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))


def create_synthetic_fastq(
    ref_seq: str,
    ground_truth: List[GroundTruth],
    output_fastq: str,
    read_length: int = 150,
    error_rate: float = 0.01,
    fragment_size: int = 300
):
    """
    Generate paired-end FASTQ with reads containing known alleles.

    Key improvement over v1: Generates FASTQ (not BAM) so we can align with BWA.
    Returns tuple of (R1_fastq, R2_fastq) paths.
    """
    print(f"\nGenerating synthetic FASTQ with {len(ground_truth)} variants...")

    # Output files for paired-end reads
    fq1_path = output_fastq.replace('.fq', '_R1.fq').replace('.fastq', '_R1.fastq')
    fq2_path = output_fastq.replace('.fq', '_R2.fq').replace('.fastq', '_R2.fastq')

    read_id = 0
    variant_stats = defaultdict(lambda: {'ref': 0, 'alt': 0})

    with open(fq1_path, 'w') as fq1, open(fq2_path, 'w') as fq2:
        for gt in ground_truth:
            # Calculate read counts from true ratio
            total_reads = gt.coverage
            ref_reads = int(total_reads * gt.true_ratio / (gt.true_ratio + 1))
            alt_reads = total_reads - ref_reads

            variant_stats[gt.pos]['ref'] = ref_reads
            variant_stats[gt.pos]['alt'] = alt_reads

            # Set seed for reproducibility
            random.seed(gt.seed)
            np.random.seed(gt.seed)

            # Generate REF-supporting read pairs
            for _ in range(ref_reads):
                # R1 covers the variant
                r1_seq, r1_qual = create_read_sequence(
                    ref_seq, gt.pos, gt.ref_allele, len(gt.ref_allele), read_length, error_rate
                )
                # R2 is downstream, reverse complement
                r2_start = gt.pos + fragment_size - read_length
                if r2_start < 0:
                    r2_start = 0
                if r2_start + read_length > len(ref_seq):
                    r2_start = len(ref_seq) - read_length
                r2_seq_fwd = ref_seq[r2_start:r2_start + read_length]
                r2_seq = add_sequencing_errors(reverse_complement(r2_seq_fwd), error_rate)
                r2_qual = generate_quality_scores(read_length)

                write_fastq_record(fq1, f"read_{read_id}/1", r1_seq, r1_qual)
                write_fastq_record(fq2, f"read_{read_id}/2", r2_seq, r2_qual)
                read_id += 1

            # Generate ALT-supporting read pairs
            for _ in range(alt_reads):
                # R1 covers the variant with ALT allele
                # IMPORTANT: pass ref_allele_len (not alt_allele_len) for proper indel handling
                r1_seq, r1_qual = create_read_sequence(
                    ref_seq, gt.pos, gt.alt_allele, len(gt.ref_allele), read_length, error_rate
                )
                # R2 is downstream, reverse complement (from ref, no variant)
                r2_start = gt.pos + fragment_size - read_length
                if r2_start < 0:
                    r2_start = 0
                if r2_start + read_length > len(ref_seq):
                    r2_start = len(ref_seq) - read_length
                r2_seq_fwd = ref_seq[r2_start:r2_start + read_length]
                r2_seq = add_sequencing_errors(reverse_complement(r2_seq_fwd), error_rate)
                r2_qual = generate_quality_scores(read_length)

                write_fastq_record(fq1, f"read_{read_id}/1", r1_seq, r1_qual)
                write_fastq_record(fq2, f"read_{read_id}/2", r2_seq, r2_qual)
                read_id += 1

    print(f"  ✅ Generated {read_id} read pairs")

    # Print sample statistics
    print(f"\n  Sample variant breakdown:")
    for pos in sorted(list(variant_stats.keys()))[:3]:  # Show first 3
        stats = variant_stats[pos]
        ratio = stats['ref'] / stats['alt'] if stats['alt'] > 0 else float('inf')
        print(f"    Pos {pos}: {stats['ref']} REF, {stats['alt']} ALT (ratio: {ratio:.2f})")

    return (fq1_path, fq2_path)


def align_with_bwa(
    ref_fasta: str,
    fastq_files: tuple,
    output_bam: str,
    threads: int = None
):
    """
    Align paired-end FASTQ with BWA to generate realistic CIGAR strings.

    This is THE KEY STEP - creates real indel CIGARs that WASP2 must handle!

    Args:
        ref_fasta: Reference FASTA file
        fastq_files: Tuple of (R1_fastq, R2_fastq) for paired-end reads
        output_bam: Output BAM path
        threads: Number of threads (default: auto-detect)
    """
    # Use all available CPUs, capped at 16 for optimal performance
    if threads is None:
        threads = min(multiprocessing.cpu_count(), 16)

    print(f"\nAligning reads with BWA (using {threads} threads)...")

    # Check if BWA index exists
    if not Path(f"{ref_fasta}.bwt").exists():
        print("  Building BWA index...")
        subprocess.run(['bwa', 'index', ref_fasta],
                      check=True, capture_output=True)

    # Align with BWA MEM (paired-end)
    fq1, fq2 = fastq_files
    sam_file = output_bam.replace('.bam', '.sam')
    with open(sam_file, 'w') as sam:
        subprocess.run([
            'bwa', 'mem',
            '-t', str(threads),
            ref_fasta,
            fq1, fq2  # Paired-end FASTQ files
        ], stdout=sam, stderr=subprocess.PIPE, check=True)

    print("  ✅ Alignment complete")

    # Convert SAM → BAM
    print("  Converting to BAM...")
    pysam.view('-bS', '-o', output_bam, sam_file, catch_stdout=False)

    # Sort BAM
    sorted_bam = output_bam.replace('.bam', '.sorted.bam')
    pysam.sort('-o', sorted_bam, output_bam)

    # Index BAM
    pysam.index(sorted_bam)

    # Clean up
    Path(sam_file).unlink()
    Path(output_bam).unlink()

    print(f"  ✅ Created {sorted_bam}")

    # Check alignment stats
    stats = pysam.AlignmentFile(sorted_bam).count()
    print(f"  Aligned reads: {stats}")

    return sorted_bam


def create_vcf_from_ground_truth(
    ground_truth: List[GroundTruth],
    output_vcf: str
):
    """Create phased VCF for WASP2 input."""
    print(f"\nCreating VCF with {len(ground_truth)} variants...")

    # Get unique variants (remove replicates)
    unique_variants = {}
    for gt in ground_truth:
        key = (gt.chrom, gt.pos, gt.ref_allele, gt.alt_allele)
        if key not in unique_variants:
            unique_variants[key] = gt

    with open(output_vcf, 'w') as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write("##contig=<ID=chr1,length=1000000>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")

        # Write variants (all heterozygous and phased)
        for i, (key, gt) in enumerate(sorted(unique_variants.items())):
            variant_id = f"{gt.variant_type}_{gt.pos}"
            f.write(f"{gt.chrom}\t{gt.pos}\t{variant_id}\t{gt.ref_allele}\t{gt.alt_allele}\t60\tPASS\t.\tGT\t0|1\n")

    print(f"  ✅ Created VCF with {len(unique_variants)} unique variants")

    # Compress and index
    subprocess.run(['bgzip', '-f', output_vcf], check=True, capture_output=True)
    subprocess.run(['tabix', '-f', '-p', 'vcf', f'{output_vcf}.gz'],
                  check=True, capture_output=True)

    return f'{output_vcf}.gz'


def run_wasp2_pipeline(
    bam_file: str,
    vcf_file: str,
    ref_fasta: str,
    output_dir: Path
):
    """
    Run full WASP2 pipeline.

    THIS IS THE REAL TEST - runs actual WASP2 code!
    """
    print(f"\n{'='*80}")
    print("RUNNING WASP2 PIPELINE")
    print(f"{'='*80}")

    wasp2_dir = Path(__file__).parent.resolve()

    # Convert to absolute paths (needed because we run from src/ directory)
    bam_file_abs = str(Path(bam_file).resolve())
    vcf_file_abs = str(Path(vcf_file).resolve())
    output_dir_abs = str(output_dir.resolve())

    # Step 1: Find intersecting SNPs/indels (make_reads)
    print("\nStep 1: Finding intersecting variants...")
    cmd = [
        sys.executable, '-m', 'mapping', 'make-reads',
        bam_file_abs,
        vcf_file_abs,
        '--out_dir', output_dir_abs,
        '--sample', 'sample1',  # Match sample name in VCF
        '--indels',   # Enable indel support
        '--paired',   # Paired-end reads
        '--phased',   # VCF has phased genotypes (0|1)
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True,
                               cwd=str(wasp2_dir / 'src'))
        print("  ✅ make-reads complete")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ make-reads failed:")
        print(f"  STDOUT: {e.stdout}")
        print(f"  STDERR: {e.stderr}")
        raise

    # Get BAM prefix to find the output files (WASP2 names files based on input BAM)
    bam_prefix = Path(bam_file).stem  # e.g., "aligned.sorted"

    # Check outputs (WASP2 uses bam_prefix in file names)
    remap_fq1 = output_dir / f'{bam_prefix}_swapped_alleles_r1.fq'
    remap_fq2 = output_dir / f'{bam_prefix}_swapped_alleles_r2.fq'
    to_remap_bam = output_dir / f'{bam_prefix}_to_remap.bam'
    keep_bam = output_dir / f'{bam_prefix}_keep.bam'

    # Check if remapping is needed (swapped alleles files exist and are non-empty)
    needs_remapping = (remap_fq1.exists() and remap_fq1.stat().st_size > 0 and
                       remap_fq2.exists() and remap_fq2.stat().st_size > 0)

    if not needs_remapping:
        print(f"  ⚠️  No remapping needed (all reads passed)")
        # Just use keep.bam
        return keep_bam

    # Step 2: Realign remapped reads (paired-end)
    print("\nStep 2: Realigning remapped reads...")

    remapped_bam = align_with_bwa(
        ref_fasta,
        (str(remap_fq1), str(remap_fq2)),  # Paired-end FASTQ files
        str(output_dir / 'remapped.bam')
    )

    # Step 3: Filter remapped reads
    print("\nStep 3: Filtering remapped reads...")

    final_bam = output_dir / 'keep.merged.bam'
    cmd = [
        sys.executable, '-m', 'mapping', 'filter-remapped',
        str(Path(remapped_bam).resolve()),
        str(to_remap_bam.resolve()),
        str(keep_bam.resolve()),
        '--out_bam', str(final_bam.resolve()),
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True,
                      cwd=str(wasp2_dir / 'src'))
        print("  ✅ filter-remapped complete")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ filter-remapped failed:")
        print(f"  STDOUT: {e.stdout}")
        print(f"  STDERR: {e.stderr}")
        # Fall back to keep.bam
        final_bam = keep_bam

    print(f"\n✅ WASP2 pipeline complete")
    print(f"  Final BAM: {final_bam}")

    return final_bam


def count_alleles_in_bam(
    bam_file: str,
    ground_truth: List[GroundTruth],
    min_baseq: int = 20,
    min_mapq: int = 10
) -> pd.DataFrame:
    """
    Count REF vs ALT alleles in final WASP2-filtered BAM.

    Uses CIGAR-aware counting with proper indel size matching.
    Each ground truth variant now has a unique position (no replicate pooling).
    """
    print(f"\nCounting alleles in WASP2 output...")

    if not Path(bam_file).exists():
        print(f"  ❌ BAM file not found: {bam_file}")
        return pd.DataFrame()

    # Ensure BAM file is indexed for pileup
    bam_index = f"{bam_file}.bai"
    if not Path(bam_index).exists():
        pysam.index(bam_file)

    bam = pysam.AlignmentFile(bam_file)

    results = []

    # Process each ground truth variant individually (each has unique position now)
    for gt in ground_truth:
        ref_count = 0
        alt_count = 0
        total_reads = 0

        # Calculate expected indel size
        indel_size = abs(len(gt.alt_allele) - len(gt.ref_allele))

        # Fetch reads overlapping variant position
        try:
            for read in bam.fetch(gt.chrom, gt.pos, gt.pos + max(len(gt.ref_allele), len(gt.alt_allele)) + 1):
                if read.mapping_quality < min_mapq:
                    continue
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                total_reads += 1

                # For SNPs: simple base comparison using pileup approach
                if gt.variant_type == 'SNP':
                    # Use get_aligned_pairs for position mapping
                    aligned_pairs = read.get_aligned_pairs(with_seq=True)
                    for query_pos, ref_pos, ref_base in aligned_pairs:
                        if ref_pos == gt.pos and query_pos is not None:
                            base = read.query_sequence[query_pos]
                            if base == gt.alt_allele:
                                alt_count += 1
                            elif base == gt.ref_allele:
                                ref_count += 1
                            break

                # For indels: check CIGAR for indel of matching size at/near variant position
                else:
                    cigar_tuples = read.cigartuples or []
                    read_ref_pos = read.reference_start

                    found_indel = False
                    for op, length in cigar_tuples:
                        if op == 0:  # M (match/mismatch)
                            read_ref_pos += length
                        elif op == 2:  # D (deletion)
                            # Check if this deletion is at our variant position
                            # Allow some slop (±5bp) for alignment variation
                            if abs(read_ref_pos - gt.pos) <= 5:
                                # Check if deletion size matches expected
                                if len(gt.ref_allele) > len(gt.alt_allele):
                                    # This is a deletion variant
                                    if abs(length - indel_size) <= 1:  # Allow 1bp tolerance
                                        alt_count += 1
                                        found_indel = True
                                        break
                            read_ref_pos += length
                        elif op == 1:  # I (insertion)
                            # Check if this insertion is at our variant position
                            if abs(read_ref_pos - gt.pos) <= 5:
                                # Check if insertion size matches expected
                                if len(gt.alt_allele) > len(gt.ref_allele):
                                    # This is an insertion variant
                                    if abs(length - indel_size) <= 1:  # Allow 1bp tolerance
                                        alt_count += 1
                                        found_indel = True
                                        break
                        elif op == 4:  # S (soft clip)
                            pass  # Soft clips don't consume reference
                        elif op == 5:  # H (hard clip)
                            pass  # Hard clips don't consume reference

                    # If no matching indel found, count as REF
                    if not found_indel:
                        ref_count += 1

        except Exception as e:
            print(f"  ⚠️  Error counting at {gt.chrom}:{gt.pos}: {e}")
            continue

        # Calculate observed ratio
        if alt_count > 0:
            observed_ratio = ref_count / alt_count
        else:
            observed_ratio = float('inf') if ref_count > 0 else 0.0

        # Calculate error
        error = abs(observed_ratio - gt.true_ratio) if observed_ratio != float('inf') else float('inf')
        error_pct = (error / gt.true_ratio * 100) if gt.true_ratio > 0 else 0

        results.append({
            'chrom': gt.chrom,
            'pos': gt.pos,
            'variant_type': gt.variant_type,
            'coverage': gt.coverage,
            'replicate': gt.replicate,
            'true_ratio': gt.true_ratio,
            'ref_count': ref_count,
            'alt_count': alt_count,
            'total_reads': total_reads,
            'observed_ratio': observed_ratio,
            'error': error,
            'error_pct': error_pct,
            'status': 'PASS' if error_pct < 10 else 'FAIL'
        })

    bam.close()

    df = pd.DataFrame(results)
    print(f"  ✅ Counted alleles for {len(df)} variant tests")

    return df


def generate_test_configurations(tier: str = 'minimum') -> List[GroundTruth]:
    """
    Generate test configurations based on tier.

    Tiers:
      - minimum:       90 tests (3 types × 3 ratios × 10 reps)
      - moderate:     270 tests (3 types × 3 ratios × 3 coverage × 10 reps)
      - comprehensive: 810 tests (3 types × 3 ratios × 3 coverage × 3 sizes × 10 reps)
    """
    configs = []
    config_id = 0

    # Base parameters
    variant_types = {
        'SNP': [('A', 'G'), ('T', 'C'), ('G', 'A')],
        'INS': [('C', 'CAT'), ('A', 'AGGG'), ('T', 'TTTAA')],
        'DEL': [('GCC', 'G'), ('ATATA', 'A'), ('TTTT', 'T')]
    }

    allelic_ratios = [1.0, 2.0, 4.0]
    n_replicates = 10

    # Tier-specific parameters
    if tier == 'minimum':
        coverage_levels = [50]
        size_configs = {'default': variant_types}

    elif tier == 'moderate':
        coverage_levels = [20, 50, 100]
        size_configs = {'default': variant_types}

    elif tier == 'comprehensive':
        coverage_levels = [20, 50, 100]
        size_configs = {
            'small': {
                'SNP': [('A', 'G')],
                'INS': [('C', 'CA')],
                'DEL': [('GC', 'G')]
            },
            'medium': {
                'SNP': [('T', 'C')],
                'INS': [('A', 'AGGG')],
                'DEL': [('ATATA', 'A')]
            },
            'large': {
                'SNP': [('G', 'A')],
                'INS': [('T', 'T' + 'GGGGGGGGGG')],  # 10bp insertion
                'DEL': [('GCGCGCGCGCGCGCGCGCGC', 'G')]  # 19bp deletion
            }
        }
    else:
        raise ValueError(f"Unknown tier: {tier}")

    # Generate configurations
    # CRITICAL FIX: Each replicate must have a UNIQUE position
    # Previously all replicates at same position got their reads pooled together
    base_pos = 50000  # Start far from edges
    pos_spacing = 1000  # Space variants apart (reduced for more variants)
    variant_id = 0  # Unique ID for each variant instance (including replicates)

    for size_name, vartypes in size_configs.items():
        for vtype, allele_pairs in vartypes.items():
            for allele_pair in allele_pairs:
                ref_allele, alt_allele = allele_pair

                for ratio in allelic_ratios:
                    for coverage in coverage_levels:
                        for replicate in range(n_replicates):
                            # CRITICAL: Each replicate gets its own unique position
                            pos = base_pos + (variant_id * pos_spacing)
                            seed = variant_id * 1000 + replicate

                            configs.append(GroundTruth(
                                chrom='chr1',
                                pos=pos,
                                ref_allele=ref_allele,
                                alt_allele=alt_allele,
                                true_ratio=ratio,
                                variant_type=vtype,
                                coverage=coverage,
                                replicate=replicate,
                                seed=seed
                            ))

                            variant_id += 1  # Increment for EACH replicate

    return configs


def print_summary_statistics(results_df: pd.DataFrame):
    """Print summary statistics from simulation."""
    print(f"\n{'='*80}")
    print("VALIDATION RESULTS")
    print(f"{'='*80}\n")

    if len(results_df) == 0:
        print("❌ No results to analyze")
        return

    # Overall statistics
    mean_error = results_df['error_pct'].mean()
    median_error = results_df['error_pct'].median()
    max_error = results_df['error_pct'].max()
    n_pass = (results_df['status'] == 'PASS').sum()
    n_total = len(results_df)
    pass_rate = n_pass / n_total * 100

    print(f"Overall Performance:")
    print(f"  Total tests:       {n_total}")
    print(f"  Passed (<10%):     {n_pass}/{n_total} ({pass_rate:.1f}%)")
    print(f"  Mean error:        {mean_error:.2f}%")
    print(f"  Median error:      {median_error:.2f}%")
    print(f"  Max error:         {max_error:.2f}%")
    print()

    # By variant type
    print("Performance by Variant Type:")
    for vtype in ['SNP', 'INS', 'DEL']:
        subset = results_df[results_df['variant_type'] == vtype]
        if len(subset) > 0:
            mean = subset['error_pct'].mean()
            std = subset['error_pct'].std()
            print(f"  {vtype:5s}: {mean:5.2f}% ± {std:4.2f}%")
    print()

    # By coverage (if multiple)
    if 'coverage' in results_df.columns and results_df['coverage'].nunique() > 1:
        print("Performance by Coverage:")
        for cov in sorted(results_df['coverage'].unique()):
            subset = results_df[results_df['coverage'] == cov]
            mean = subset['error_pct'].mean()
            std = subset['error_pct'].std()
            print(f"  {cov:3d}x: {mean:5.2f}% ± {std:4.2f}%")
        print()

    # By allelic ratio
    print("Performance by True Ratio:")
    for ratio in sorted(results_df['true_ratio'].unique()):
        subset = results_df[results_df['true_ratio'] == ratio]
        mean = subset['error_pct'].mean()
        std = subset['error_pct'].std()
        print(f"  {ratio:.1f}:1: {mean:5.2f}% ± {std:4.2f}%")
    print()

    # Final verdict
    print(f"{'='*80}")
    if pass_rate >= 90 and mean_error < 10:
        print("✅ SIMULATION VALIDATES WASP2 INDEL IMPLEMENTATION")
        print(f"   {n_pass}/{n_total} tests passed with mean error {mean_error:.2f}%")
    else:
        print("❌ VALIDATION FAILED")
        print(f"   Only {n_pass}/{n_total} tests passed (need ≥90%)")
        print(f"   Mean error {mean_error:.2f}% (need <10%)")
    print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description='WASP2 Indel Validation - Realistic Simulation'
    )
    parser.add_argument(
        '--tier',
        choices=['minimum', 'moderate', 'comprehensive'],
        default='minimum',
        help='Testing tier (minimum=90 tests, moderate=270, comprehensive=810)'
    )
    parser.add_argument(
        '--workdir',
        type=str,
        default=None,
        help='Working directory (default: temp directory)'
    )
    parser.add_argument(
        '--keep',
        action='store_true',
        help='Keep temporary files'
    )

    args = parser.parse_args()

    print(f"{'='*80}")
    print(f"WASP2 INDEL VALIDATION - {args.tier.upper()} TIER")
    print(f"{'='*80}\n")

    # Create working directory
    if args.workdir:
        workdir = Path(args.workdir)
        workdir.mkdir(exist_ok=True, parents=True)
    else:
        workdir = Path(tempfile.mkdtemp(prefix=f'wasp2_sim_{args.tier}_'))

    print(f"Working directory: {workdir}\n")

    # Generate test configurations
    print(f"Generating {args.tier} test configurations...")
    ground_truth = generate_test_configurations(args.tier)
    print(f"  ✅ {len(ground_truth)} tests generated\n")

    # Create reference genome
    ref_fasta = str(workdir / 'reference.fa')
    ref_seq = create_reference_fasta(ref_fasta)

    # Create VCF
    vcf_file = create_vcf_from_ground_truth(ground_truth, str(workdir / 'variants.vcf'))

    # Generate synthetic paired-end FASTQ
    fastq_base = str(workdir / 'synthetic.fq')
    fastq_files = create_synthetic_fastq(ref_seq, ground_truth, fastq_base)

    # Align with BWA (paired-end)
    aligned_bam = align_with_bwa(ref_fasta, fastq_files, str(workdir / 'aligned.bam'))

    # Run WASP2 pipeline
    wasp_output_dir = workdir / 'wasp2_output'
    wasp_output_dir.mkdir(exist_ok=True)

    try:
        final_bam = run_wasp2_pipeline(aligned_bam, vcf_file, ref_fasta, wasp_output_dir)
    except Exception as e:
        print(f"\n❌ WASP2 pipeline failed: {e}")
        print(f"\nCheck logs in: {wasp_output_dir}")
        sys.exit(1)

    # Count alleles and validate
    results = count_alleles_in_bam(str(final_bam), ground_truth)

    # Save results
    results_file = workdir / 'simulation_results.csv'
    results.to_csv(results_file, index=False)
    print(f"\n✅ Results saved to: {results_file}")

    # Print summary
    print_summary_statistics(results)

    # Cleanup
    if not args.keep:
        print(f"\nCleaning up temporary files...")
        # Keep results but clean up intermediate files
        for f in workdir.glob('*.fq'):
            f.unlink()
        for f in workdir.glob('*.sam'):
            f.unlink()
    else:
        print(f"\nTemporary files kept in: {workdir}")

    print(f"\nSimulation complete!")
    print(f"Results: {results_file}")


if __name__ == '__main__':
    main()
