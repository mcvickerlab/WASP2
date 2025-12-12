#!/usr/bin/env python3
"""
Detailed analysis of WASP2-Rust vs Python DEV divergence.
Traces specific reads through both pipelines to identify root cause.
"""

import sys
import json
import pysam
import pandas as pd
from pathlib import Path
from collections import defaultdict

WORKDIR = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/atac_gm12878/chr22_comparison")
RUST_DIR = WORKDIR / "rust"
PYTHON_DIR = WORKDIR / "python"
ANALYSIS_DIR = WORKDIR / "analysis"

def parse_wasp_name(name):
    """Parse WASP-encoded read name.
    Format: {original_name}_WASP_{r1_pos}_{r2_pos}_{idx}_{total}
    """
    parts = name.split("_WASP_")
    if len(parts) != 2:
        return None

    original = parts[0]
    wasp_parts = parts[1].split("_")
    if len(wasp_parts) < 4:
        return None

    return {
        "original_name": original,
        "r1_pos": int(wasp_parts[0]),
        "r2_pos": int(wasp_parts[1]),
        "hap_idx": int(wasp_parts[2]),
        "hap_total": int(wasp_parts[3])
    }


def analyze_step1_bed_files():
    """Compare VCF to BED conversion outputs."""
    print("\n" + "="*60)
    print("STEP 1 ANALYSIS: VCF to BED")
    print("="*60)

    rust_bed = RUST_DIR / "step1_vcf_to_bed" / "variants.bed"
    python_bed = PYTHON_DIR / "step1_vcf_to_bed" / "variants.bed"

    if not rust_bed.exists() or not python_bed.exists():
        print("BED files not found. Run pipelines first.")
        return

    # Read both BED files
    rust_variants = set()
    with open(rust_bed) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                rust_variants.add((parts[0], parts[1], parts[2]))

    python_variants = set()
    with open(python_bed) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                python_variants.add((parts[0], parts[1], parts[2]))

    print(f"Rust variants:   {len(rust_variants)}")
    print(f"Python variants: {len(python_variants)}")
    print(f"Common:          {len(rust_variants & python_variants)}")
    print(f"Rust-only:       {len(rust_variants - python_variants)}")
    print(f"Python-only:     {len(python_variants - rust_variants)}")

    # Key difference: Python uses '-v snps' to filter to SNPs only
    print("\nNOTE: Python DEV uses 'bcftools view -v snps' to filter to SNPs only!")
    print("Rust may include indels depending on vcf_to_bed implementation.")


def analyze_step3_make_reads():
    """Compare make-reads outputs - the critical step."""
    print("\n" + "="*60)
    print("STEP 3 ANALYSIS: Make Reads (CRITICAL)")
    print("="*60)

    rust_wasp = RUST_DIR / "step3_make_reads" / "wasp_names.txt"
    python_wasp = PYTHON_DIR / "step3_make_reads" / "wasp_names.txt"

    if not rust_wasp.exists() or not python_wasp.exists():
        print("WASP name files not found. Run pipelines first.")
        return

    # Parse WASP names
    rust_reads = {}
    with open(rust_wasp) as f:
        for line in f:
            name = line.strip()
            parsed = parse_wasp_name(name)
            if parsed:
                orig = parsed["original_name"]
                if orig not in rust_reads:
                    rust_reads[orig] = []
                rust_reads[orig].append(parsed)

    python_reads = {}
    with open(python_wasp) as f:
        for line in f:
            name = line.strip()
            parsed = parse_wasp_name(name)
            if parsed:
                orig = parsed["original_name"]
                if orig not in python_reads:
                    python_reads[orig] = []
                python_reads[orig].append(parsed)

    print(f"Rust unique original reads:   {len(rust_reads)}")
    print(f"Python unique original reads: {len(python_reads)}")

    # Compare haplotype counts per read
    rust_hap_counts = defaultdict(int)
    python_hap_counts = defaultdict(int)

    for orig, haps in rust_reads.items():
        rust_hap_counts[len(haps)] += 1
    for orig, haps in python_reads.items():
        python_hap_counts[len(haps)] += 1

    print("\nHaplotypes per read distribution:")
    print("Count | Rust   | Python")
    all_counts = sorted(set(rust_hap_counts.keys()) | set(python_hap_counts.keys()))
    for c in all_counts:
        print(f"  {c}   | {rust_hap_counts.get(c, 0):6d} | {python_hap_counts.get(c, 0):6d}")

    # Find reads with different haplotype counts
    common_reads = set(rust_reads.keys()) & set(python_reads.keys())
    diff_hap_count = []
    for orig in common_reads:
        rust_count = len(rust_reads[orig])
        python_count = len(python_reads[orig])
        if rust_count != python_count:
            diff_hap_count.append({
                "read": orig,
                "rust_haps": rust_count,
                "python_haps": python_count
            })

    print(f"\nReads with different haplotype counts: {len(diff_hap_count)}")
    if diff_hap_count:
        print("First 10 examples:")
        for d in diff_hap_count[:10]:
            print(f"  {d['read']}: Rust={d['rust_haps']}, Python={d['python_haps']}")

        # Save for further analysis
        pd.DataFrame(diff_hap_count).to_csv(ANALYSIS_DIR / "diff_haplotype_counts.csv", index=False)


def analyze_step5_filter():
    """Compare filter-remapped outputs."""
    print("\n" + "="*60)
    print("STEP 5 ANALYSIS: Filter Remapped (CRITICAL)")
    print("="*60)

    rust_kept = ANALYSIS_DIR / "rust_kept_sorted.txt"
    python_kept = ANALYSIS_DIR / "python_kept_sorted.txt"
    python_keeps_rust_removes = ANALYSIS_DIR / "step5_python_keeps_rust_removes.txt"

    if not python_keeps_rust_removes.exists():
        print("Comparison files not found. Run 04_compare_outputs.sh first.")
        return

    # Load reads that Python keeps but Rust removes
    divergent_reads = []
    with open(python_keeps_rust_removes) as f:
        divergent_reads = [line.strip() for line in f if line.strip()]

    print(f"Reads Python keeps but Rust removes: {len(divergent_reads)}")

    if not divergent_reads:
        print("No divergent reads found!")
        return

    # Sample some divergent reads for detailed analysis
    sample_reads = divergent_reads[:100]

    print(f"\nAnalyzing {len(sample_reads)} sample divergent reads...")

    # Analyze in remapped BAM files
    rust_remapped = RUST_DIR / "step4_remap" / "remapped.bam"
    python_remapped = PYTHON_DIR / "step4_remap" / "remapped.bam"

    if rust_remapped.exists() and python_remapped.exists():
        analysis_results = []

        # Create index of sample reads for quick lookup
        sample_set = set(sample_reads)

        # Check Rust remapped
        with pysam.AlignmentFile(str(rust_remapped), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                # Parse WASP name to get original name
                parsed = parse_wasp_name(read.query_name)
                if parsed and parsed["original_name"] in sample_set:
                    analysis_results.append({
                        "original_name": parsed["original_name"],
                        "wasp_name": read.query_name,
                        "source": "rust",
                        "chrom": read.reference_name,
                        "pos": read.reference_start,
                        "mate_pos": read.next_reference_start,
                        "wasp_r1_pos": parsed["r1_pos"],
                        "wasp_r2_pos": parsed["r2_pos"],
                        "pos_match_r1": read.reference_start == parsed["r1_pos"],
                        "pos_match_r2": read.next_reference_start == parsed["r2_pos"],
                    })

        # Check Python remapped
        with pysam.AlignmentFile(str(python_remapped), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                parsed = parse_wasp_name(read.query_name)
                if parsed and parsed["original_name"] in sample_set:
                    analysis_results.append({
                        "original_name": parsed["original_name"],
                        "wasp_name": read.query_name,
                        "source": "python",
                        "chrom": read.reference_name,
                        "pos": read.reference_start,
                        "mate_pos": read.next_reference_start,
                        "wasp_r1_pos": parsed["r1_pos"],
                        "wasp_r2_pos": parsed["r2_pos"],
                        "pos_match_r1": read.reference_start == parsed["r1_pos"],
                        "pos_match_r2": read.next_reference_start == parsed["r2_pos"],
                    })

        if analysis_results:
            df = pd.DataFrame(analysis_results)
            df.to_csv(ANALYSIS_DIR / "divergent_reads_analysis.csv", index=False)

            # Summarize position match patterns
            print("\nPosition match summary for divergent reads:")
            summary = df.groupby(["source", "pos_match_r1", "pos_match_r2"]).size()
            print(summary)


def main():
    """Run all analyses."""
    print("WASP2-Rust vs Python DEV Detailed Analysis")
    print("="*60)

    ANALYSIS_DIR.mkdir(exist_ok=True)

    analyze_step1_bed_files()
    analyze_step3_make_reads()
    analyze_step5_filter()

    print("\n" + "="*60)
    print("Analysis complete!")
    print(f"Results saved to: {ANALYSIS_DIR}")
    print("="*60)


if __name__ == "__main__":
    main()
