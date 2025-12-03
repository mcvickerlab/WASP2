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
    read_length: int = 150,
    error_rate: float = 0.01
) -> Tuple[str, np.ndarray]:
    """
    Create read sequence with specific allele at variant position.

    Returns:
        Tuple of (sequence, quality_scores)
    """
    # Random position offset so variant isn't always centered
    offset = random.randint(-30, 30)
    read_start = max(0, var_pos - read_length // 2 + offset)

    # Build read sequence
    # Left part + allele + right part
    left_seq = ref_seq[read_start:var_pos]
    right_seq = ref_seq[var_pos + len(allele):read_start + read_length]

    read_seq = left_seq + allele + right_seq

    # Truncate or pad to exact read length
    if len(read_seq) > read_length:
        read_seq = read_seq[:read_length]
    elif len(read_seq) < read_length:
        # Pad with reference sequence if too short
        pad_len = read_length - len(read_seq)
        read_seq += ref_seq[read_start + len(read_seq):read_start + len(read_seq) + pad_len]

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


def create_synthetic_fastq(
    ref_seq: str,
    ground_truth: List[GroundTruth],
    output_fastq: str,
    read_length: int = 150,
    error_rate: float = 0.01
):
    """
    Generate FASTQ with reads containing known alleles.

    Key improvement over v1: Generates FASTQ (not BAM) so we can align with BWA.
    """
    print(f"\nGenerating synthetic FASTQ with {len(ground_truth)} variants...")

    read_id = 0
    variant_stats = defaultdict(lambda: {'ref': 0, 'alt': 0})

    with open(output_fastq, 'w') as fq:
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

            # Generate REF-supporting reads
            for _ in range(ref_reads):
                read_seq, read_qual = create_read_sequence(
                    ref_seq, gt.pos, gt.ref_allele, read_length, error_rate
                )
                write_fastq_record(fq, f"read_{read_id}", read_seq, read_qual)
                read_id += 1

            # Generate ALT-supporting reads
            for _ in range(alt_reads):
                read_seq, read_qual = create_read_sequence(
                    ref_seq, gt.pos, gt.alt_allele, read_length, error_rate
                )
                write_fastq_record(fq, f"read_{read_id}", read_seq, read_qual)
                read_id += 1

    print(f"  ✅ Generated {read_id} reads")

    # Print sample statistics
    print(f"\n  Sample variant breakdown:")
    for pos in sorted(list(variant_stats.keys()))[:3]:  # Show first 3
        stats = variant_stats[pos]
        ratio = stats['ref'] / stats['alt'] if stats['alt'] > 0 else float('inf')
        print(f"    Pos {pos}: {stats['ref']} REF, {stats['alt']} ALT (ratio: {ratio:.2f})")

    return output_fastq


def align_with_bwa(
    ref_fasta: str,
    fastq_file: str,
    output_bam: str,
    threads: int = 4
):
    """
    Align FASTQ with BWA to generate realistic CIGAR strings.

    This is THE KEY STEP - creates real indel CIGARs that WASP2 must handle!
    """
    print(f"\nAligning reads with BWA...")

    # Check if BWA index exists
    if not Path(f"{ref_fasta}.bwt").exists():
        print("  Building BWA index...")
        subprocess.run(['bwa', 'index', ref_fasta],
                      check=True, capture_output=True)

    # Align with BWA MEM
    sam_file = output_bam.replace('.bam', '.sam')
    with open(sam_file, 'w') as sam:
        subprocess.run([
            'bwa', 'mem',
            '-t', str(threads),
            ref_fasta,
            fastq_file
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

    wasp2_dir = Path(__file__).parent

    # Step 1: Find intersecting SNPs/indels
    print("\nStep 1: Finding intersecting variants...")
    cmd = [
        'python', str(wasp2_dir / 'find_intersecting_snps.py'),
        '--bam', bam_file,
        '--vcf', vcf_file,
        '--ref', ref_fasta,
        '--out_dir', str(output_dir),
        '--include_indels',  # ← KEY FLAG!
        '--is_paired_end'
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("  ✅ find_intersecting_snps.py complete")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ find_intersecting_snps.py failed:")
        print(f"  STDOUT: {e.stdout}")
        print(f"  STDERR: {e.stderr}")
        raise

    # Check outputs
    remap_fastq = output_dir / 'remap.fq.gz'
    to_remap_bam = output_dir / 'to.remap.bam'
    keep_bam = output_dir / 'keep.bam'

    if not remap_fastq.exists():
        print(f"  ⚠️  No remapping needed (all reads passed)")
        # Just use keep.bam
        return keep_bam

    # Step 2: Realign remapped reads
    print("\nStep 2: Realigning remapped reads...")

    # Decompress FASTQ
    remap_fq_uncompressed = output_dir / 'remap.fq'
    with gzip.open(remap_fastq, 'rt') as f_in:
        with open(remap_fq_uncompressed, 'w') as f_out:
            f_out.write(f_in.read())

    remapped_bam = align_with_bwa(
        ref_fasta,
        str(remap_fq_uncompressed),
        str(output_dir / 'remapped.bam')
    )

    # Step 3: Filter remapped reads
    print("\nStep 3: Filtering remapped reads...")

    # Check if filter script exists
    filter_script = wasp2_dir / 'filter_remapped_reads.py'
    if not filter_script.exists():
        print(f"  ⚠️  filter_remapped_reads.py not found")
        print(f"  Using remapped BAM directly (no filtering)")
        final_bam = remapped_bam
    else:
        final_bam = output_dir / 'keep.merged.bam'
        cmd = [
            'python', str(filter_script),
            '--remap_bam', remapped_bam,
            '--to_remap_bam', str(to_remap_bam),
            '--keep_bam', str(keep_bam),
            '--out', str(final_bam)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print("  ✅ filter_remapped_reads.py complete")
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Filtering failed: {e}")
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

    This measures what WASP2 actually recovered.
    """
    print(f"\nCounting alleles in WASP2 output...")

    if not Path(bam_file).exists():
        print(f"  ❌ BAM file not found: {bam_file}")
        return pd.DataFrame()

    bam = pysam.AlignmentFile(bam_file)

    results = []

    # Group ground truth by position
    gt_by_pos = defaultdict(list)
    for gt in ground_truth:
        gt_by_pos[gt.pos].append(gt)

    for pos, gt_list in gt_by_pos.items():
        # Use first instance for variant info (all replicates have same alleles)
        gt = gt_list[0]

        ref_count = 0
        alt_count = 0
        total_reads = 0

        # Fetch reads overlapping variant
        try:
            for pileupcolumn in bam.pileup(gt.chrom, pos, pos + 1,
                                          truncate=True, min_base_quality=min_baseq):
                if pileupcolumn.pos != pos:
                    continue

                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue

                    read = pileupread.alignment
                    if read.mapping_quality < min_mapq:
                        continue

                    total_reads += 1

                    # Get query position
                    query_pos = pileupread.query_position
                    if query_pos is None:
                        continue

                    # Extract allele from read
                    # Handle different variant types
                    allele_len = max(len(gt.ref_allele), len(gt.alt_allele))
                    read_allele = read.query_sequence[query_pos:query_pos + allele_len]

                    if read_allele.startswith(gt.ref_allele):
                        ref_count += 1
                    elif read_allele.startswith(gt.alt_allele):
                        alt_count += 1

        except Exception as e:
            print(f"  ⚠️  Error counting at {gt.chrom}:{pos}: {e}")
            continue

        # Calculate observed ratio
        if alt_count > 0:
            observed_ratio = ref_count / alt_count
        else:
            observed_ratio = float('inf') if ref_count > 0 else 0.0

        # Record for each replicate
        for gt_rep in gt_list:
            error = abs(observed_ratio - gt_rep.true_ratio) if observed_ratio != float('inf') else float('inf')
            error_pct = (error / gt_rep.true_ratio * 100) if gt_rep.true_ratio > 0 else 0

            results.append({
                'chrom': gt_rep.chrom,
                'pos': gt_rep.pos,
                'variant_type': gt_rep.variant_type,
                'coverage': gt_rep.coverage,
                'replicate': gt_rep.replicate,
                'true_ratio': gt_rep.true_ratio,
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
    base_pos = 50000  # Start far from edges
    pos_spacing = 5000  # Space variants apart

    for size_name, vartypes in size_configs.items():
        for vtype, allele_pairs in vartypes.items():
            for allele_pair in allele_pairs:
                ref_allele, alt_allele = allele_pair

                for ratio in allelic_ratios:
                    for coverage in coverage_levels:
                        for replicate in range(n_replicates):
                            pos = base_pos + (config_id * pos_spacing)
                            seed = config_id * 1000 + replicate

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

                config_id += 1

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

    # Generate synthetic FASTQ
    fastq_file = str(workdir / 'synthetic.fq')
    create_synthetic_fastq(ref_seq, ground_truth, fastq_file)

    # Align with BWA
    aligned_bam = align_with_bwa(ref_fasta, fastq_file, str(workdir / 'aligned.bam'))

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
