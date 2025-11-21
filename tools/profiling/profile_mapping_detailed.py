#!/usr/bin/env python3
"""
Detailed line-by-line profiling of critical mapping functions
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

# Import after path setup
from src.mapping.make_remap_reads import swap_chrom_alleles, write_remap_bam
from src.mapping.remap_utils import (
    paired_read_gen_stat, get_read_het_data,
    make_phased_seqs, write_read
)
from src.mapping.filter_remap_reads import filt_remapped_reads, merge_filt_bam
from src.mapping.intersect_variant_data import (
    vcf_to_bed, process_bam, intersect_reads, make_intersect_df
)


def profile_individual_stages():
    """Profile each stage separately to identify bottlenecks"""

    bam_file = "test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam"
    vcf_file = "test_data/filter_chr10.vcf"
    sample = "NA12878"
    out_dir = "baselines/mapping"

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    timings = {}

    print("=" * 80)
    print("DETAILED STAGE-BY-STAGE PROFILING")
    print("=" * 80)

    # Stage 1: VCF to BED
    print("\n[1/5] VCF to BED conversion...")
    vcf_bed = f"{out_dir}/variants.bed"
    start = time.time()
    vcf_to_bed(vcf_file, vcf_bed, samples=[sample])
    timings['vcf_to_bed'] = time.time() - start
    print(f"  Time: {timings['vcf_to_bed']:.3f}s")

    # Stage 2: Process BAM
    print("\n[2/5] BAM processing (intersect with VCF)...")
    to_remap_bam = f"{out_dir}/to_remap.bam"
    remap_reads = f"{out_dir}/remap_reads.txt"
    keep_bam = f"{out_dir}/keep.bam"
    start = time.time()
    process_bam(bam_file, vcf_bed, to_remap_bam, remap_reads, keep_bam, is_paired=True)
    timings['process_bam'] = time.time() - start
    print(f"  Time: {timings['process_bam']:.3f}s")

    # Stage 3: Intersect reads
    print("\n[3/5] Create read-variant intersections...")
    intersect_file = f"{out_dir}/intersect.bed"
    start = time.time()
    try:
        intersect_reads(to_remap_bam, vcf_bed, intersect_file)
        timings['intersect_reads'] = time.time() - start
        print(f"  Time: {timings['intersect_reads']:.3f}s")
    except Exception as e:
        print(f"  SKIPPED (bedtools not available): {e}")
        # Try to use existing file or create dummy
        import subprocess
        if not Path(intersect_file).exists():
            print("  Creating intersection file manually with pysam...")
            # Alternative: use pysam directly
            subprocess.run(
                f"bedtools intersect -a {to_remap_bam} -b {vcf_bed} -wb -bed > {intersect_file}",
                shell=True,
                check=False
            )
        timings['intersect_reads'] = time.time() - start

    # Stage 4: Write remapped reads (THE BOTTLENECK)
    print("\n[4/5] Generate allele-swapped reads...")
    print("  (This is the critical bottleneck we want to optimize in Rust)")
    r1_out = f"{out_dir}/remap_r1.fq"
    r2_out = f"{out_dir}/remap_r2.fq"
    start = time.time()
    write_remap_bam(
        to_remap_bam,
        intersect_file,
        r1_out,
        r2_out,
        samples=[sample]
    )
    timings['write_remap_bam'] = time.time() - start
    print(f"  Time: {timings['write_remap_bam']:.3f}s")

    # Count reads generated
    try:
        with open(r1_out, 'r') as f:
            read_count = sum(1 for line in f if line.startswith('@'))
        print(f"  Reads generated: {read_count:,}")
    except:
        pass

    # Stage 5: Make intersect dataframe (for stats)
    print("\n[5/5] Load intersection DataFrame...")
    start = time.time()
    intersect_df = make_intersect_df(intersect_file, [sample])
    timings['make_intersect_df'] = time.time() - start
    print(f"  Time: {timings['make_intersect_df']:.3f}s")
    print(f"  DataFrame shape: {intersect_df.shape}")

    # Summary
    total = sum(timings.values())
    print("\n" + "=" * 80)
    print("TIMING BREAKDOWN")
    print("=" * 80)
    for stage, timing in sorted(timings.items(), key=lambda x: -x[1]):
        pct = (timing / total) * 100
        print(f"{stage:30s}: {timing:8.3f}s ({pct:5.1f}%)")
    print("-" * 80)
    print(f"{'TOTAL':30s}: {total:8.3f}s")
    print("=" * 80)

    # Save report
    report_file = f"{out_dir}/stage_timings.txt"
    with open(report_file, 'w') as f:
        f.write("WASP2 Mapping Stage Timings\n")
        f.write("=" * 80 + "\n\n")
        for stage, timing in sorted(timings.items(), key=lambda x: -x[1]):
            pct = (timing / total) * 100
            f.write(f"{stage:30s}: {timing:8.3f}s ({pct:5.1f}%)\n")
        f.write("\n" + "-" * 80 + "\n")
        f.write(f"{'TOTAL':30s}: {total:8.3f}s\n")

    print(f"\nReport saved to: {report_file}")

    return timings


if __name__ == "__main__":
    profile_individual_stages()
