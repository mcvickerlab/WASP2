#!/usr/bin/env python3
"""
Analyze SNV vs INDEL read counts before/after WASP filtering.

For Figure 1C: Shows the proportion of reads overlapping SNVs vs INDELs
that pass the WASP remapping filter.

Usage:
    python analyze_snv_vs_indel.py <results_dir>
"""

import sys
import subprocess
import json
from pathlib import Path
import tempfile
import os

# Use conda-installed tools
CONDA_BIN = "/iblm/netapp/home/jjaureguy/mambaforge/envs/WASP2_dev2/bin"
BEDTOOLS = f"{CONDA_BIN}/bedtools"
SAMTOOLS = f"{CONDA_BIN}/samtools"


def parse_bed_by_variant_type(bed_path):
    """
    Parse BED file and separate SNVs from INDELs.

    BED format: chr, start, end, ref, alt, genotype
    SNV: ref and alt are single characters of same length
    INDEL: ref and alt have different lengths
    """
    snvs = []
    indels = []

    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom, start, end, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]

            # SNV if ref and alt are both single bases
            if len(ref) == 1 and len(alt) == 1:
                snvs.append((chrom, start, end, ref, alt))
            else:
                # INDEL (insertion or deletion)
                indels.append((chrom, start, end, ref, alt))

    return snvs, indels


def write_bed(variants, output_path):
    """Write variants to BED file."""
    with open(output_path, 'w') as f:
        for chrom, start, end, ref, alt in variants:
            f.write(f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\n")


def count_reads_overlapping_bed(bam_path, bed_path):
    """Count unique reads overlapping positions in BED file."""
    # Use bedtools intersect to get reads overlapping variants
    cmd = f"{BEDTOOLS} intersect -a {bam_path} -b {bed_path} -u | {SAMTOOLS} view -c"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip()) if result.stdout.strip() else 0


def count_reads_overlapping_bed_pairs(bam_path, bed_path):
    """Count unique read pairs (fragments) overlapping positions in BED file."""
    # Use bedtools intersect and count unique read names
    cmd = f"{BEDTOOLS} intersect -a {bam_path} -b {bed_path} -u | {SAMTOOLS} view | cut -f1 | sort -u | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip()) if result.stdout.strip() else 0


def analyze_results_dir(results_dir):
    """Analyze a WASP2-Rust results directory for SNV vs INDEL stats."""
    results_dir = Path(results_dir)

    # Find the BED file
    bed_files = list(results_dir.glob("*.bed"))
    if not bed_files:
        print(f"No BED file found in {results_dir}")
        return None
    bed_path = bed_files[0]

    # Find the BAM files
    original_bam = list(results_dir.glob("*_sorted.bam")) or list(results_dir.glob("subset_*.bam"))
    keep_bam = results_dir / "remap_keep.bam"
    remapped_bam = results_dir / "remapped.bam"

    if not original_bam:
        # Try to find input BAM from stats
        stats_file = results_dir / "unified_stats.json"
        if stats_file.exists():
            with open(stats_file) as f:
                stats = json.load(f)
                original_bam = [Path(stats.get('bam_file', ''))]

    if not original_bam or not original_bam[0].exists():
        print(f"Original BAM not found in {results_dir}")
        return None

    original_bam = original_bam[0]
    print(f"Using BAM: {original_bam}")
    print(f"Using BED: {bed_path}")

    # Parse BED to separate SNVs and INDELs
    snvs, indels = parse_bed_by_variant_type(bed_path)
    print(f"Found {len(snvs)} SNVs and {len(indels)} INDELs in BED file")

    if len(indels) == 0:
        print("No INDELs found - this run may not include INDELs")
        return None

    # Create temporary BED files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        snv_bed = f.name
        write_bed(snvs, snv_bed)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        indel_bed = f.name
        write_bed(indels, indel_bed)

    print("\nCounting reads overlapping variants...")

    # Count reads in original BAM overlapping each type
    snv_reads_original = count_reads_overlapping_bed_pairs(original_bam, snv_bed)
    indel_reads_original = count_reads_overlapping_bed_pairs(original_bam, indel_bed)

    print(f"Original BAM - SNV overlapping read pairs: {snv_reads_original}")
    print(f"Original BAM - INDEL overlapping read pairs: {indel_reads_original}")

    # Count reads in keep BAM overlapping each type
    if keep_bam.exists():
        snv_reads_kept = count_reads_overlapping_bed_pairs(keep_bam, snv_bed)
        indel_reads_kept = count_reads_overlapping_bed_pairs(keep_bam, indel_bed)

        print(f"Keep BAM - SNV overlapping read pairs: {snv_reads_kept}")
        print(f"Keep BAM - INDEL overlapping read pairs: {indel_reads_kept}")

        # Calculate pass rates
        snv_pass_rate = snv_reads_kept / snv_reads_original * 100 if snv_reads_original > 0 else 0
        indel_pass_rate = indel_reads_kept / indel_reads_original * 100 if indel_reads_original > 0 else 0

        print(f"\n=== WASP Filter Pass Rates ===")
        print(f"SNV reads pass rate: {snv_pass_rate:.1f}%")
        print(f"INDEL reads pass rate: {indel_pass_rate:.1f}%")

        results = {
            'snv_count': len(snvs),
            'indel_count': len(indels),
            'snv_reads_original': snv_reads_original,
            'indel_reads_original': indel_reads_original,
            'snv_reads_kept': snv_reads_kept,
            'indel_reads_kept': indel_reads_kept,
            'snv_pass_rate': snv_pass_rate,
            'indel_pass_rate': indel_pass_rate,
        }

        # Save results
        output_file = results_dir / "snv_vs_indel_analysis.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_file}")

        return results
    else:
        print(f"Keep BAM not found: {keep_bam}")
        return None


if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Default to most recent indel results
        results_dir = "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/results/unified_15m_indels_2025-12-04_04-06-10"
    else:
        results_dir = sys.argv[1]

    analyze_results_dir(results_dir)
