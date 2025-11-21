#!/usr/bin/env python3
"""
Profile the critical allele-swapping bottleneck
This is the function we want to rewrite in Rust
"""
import sys
import time
import cProfile
import pstats
from pathlib import Path
from pstats import SortKey

sys.path.insert(0, str(Path(__file__).parent / "src"))

# Create a minimal intersection file for testing
def create_test_intersect():
    """Create intersection file using pysam directly (no bedtools needed)"""
    import pysam

    bam_file = "baselines/mapping/to_remap.bam"
    vcf_bed = "baselines/mapping/variants.bed"
    out_file = "baselines/mapping/intersect.bed"

    print("Creating intersection file with pysam...")

    # Load VCF bed
    variants = {}
    with open(vcf_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chrom not in variants:
                variants[chrom] = []
            variants[chrom].append((start, end, line.strip()))

    # Intersect with BAM
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(out_file, "w") as out:
        for read in bam.fetch():
            if read.is_unmapped:
                continue

            chrom = read.reference_name
            read_start = read.reference_start
            read_end = read.reference_end

            # Check overlap with variants
            if chrom in variants:
                for v_start, v_end, v_line in variants[chrom]:
                    if read_start < v_end and read_end > v_start:
                        # Overlap found
                        mate = "1" if read.is_read1 else "2"
                        out.write(f"{chrom}\t{read_start}\t{read_end}\t"
                                f"{read.query_name}/{mate}\t{read.mapping_quality}\t"
                                f"{'+'if not read.is_reverse else '-'}\t"
                                f"{v_line}\n")

    print(f"  Intersection file created: {out_file}")
    return out_file


def profile_write_remap():
    """Profile the write_remap_bam function - the critical bottleneck"""
    from src.mapping.make_remap_reads import write_remap_bam

    # Input files (already created from previous run)
    to_remap_bam = "baselines/mapping/to_remap.bam"
    intersect_file = "baselines/mapping/intersect.bed"

    # Check if intersect file exists, if not create it
    if not Path(intersect_file).exists():
        print("Intersection file missing, creating it...")
        intersect_file = create_test_intersect()

    # Output files
    r1_out = "baselines/mapping/remap_r1.fq"
    r2_out = "baselines/mapping/remap_r2.fq"
    sample = ["NA12878"]

    print("=" * 80)
    print("PROFILING ALLELE SWAPPING (write_remap_bam)")
    print("=" * 80)
    print(f"Input BAM: {to_remap_bam}")
    print(f"Intersect: {intersect_file}")
    print(f"Sample: {sample}")
    print()

    # Profile with cProfile
    profiler = cProfile.Profile()

    print("Starting detailed profiling...")
    profiler.enable()

    start = time.time()
    try:
        write_remap_bam(
            to_remap_bam,
            intersect_file,
            r1_out,
            r2_out,
            sample,
            max_seqs=64
        )
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
    finally:
        profiler.disable()

    elapsed = time.time() - start

    print()
    print("=" * 80)
    print("TIMING RESULTS")
    print("=" * 80)
    print(f"Total time for allele swapping: {elapsed:.3f}s")
    print()

    # Save detailed profile
    profile_file = "baselines/mapping/allele_swap_profile.txt"
    with open(profile_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("ALLELE SWAPPING PROFILE (write_remap_bam)\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total time: {elapsed:.3f}s\n\n")

        # Top functions by cumulative time
        f.write("TOP 50 FUNCTIONS BY CUMULATIVE TIME\n")
        f.write("=" * 80 + "\n")
        stats = pstats.Stats(profiler, stream=f)
        stats.strip_dirs()
        stats.sort_stats(SortKey.CUMULATIVE)
        stats.print_stats(50)

        # Top functions by self time
        f.write("\n" + "=" * 80 + "\n")
        f.write("TOP 50 FUNCTIONS BY SELF TIME\n")
        f.write("=" * 80 + "\n")
        stats.sort_stats(SortKey.TIME)
        stats.print_stats(50)

        # Most called functions
        f.write("\n" + "=" * 80 + "\n")
        f.write("TOP 30 MOST CALLED FUNCTIONS\n")
        f.write("=" * 80 + "\n")
        stats.sort_stats(SortKey.CALLS)
        stats.print_stats(30)

    print(f"Detailed profile saved to: {profile_file}")

    # Print summary to console
    print("\n" + "=" * 80)
    print("TOP 20 HOTSPOTS (by cumulative time)")
    print("=" * 80)
    stats = pstats.Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats(SortKey.CUMULATIVE)
    stats.print_stats(20)

    print("\n" + "=" * 80)
    print("TOP 20 HOTSPOTS (by self time)")
    print("=" * 80)
    stats.sort_stats(SortKey.TIME)
    stats.print_stats(20)

    # Count output
    try:
        with open(r1_out) as f:
            reads = sum(1 for line in f if line.startswith('@'))
        print(f"\nReads generated: {reads:,}")
        print(f"Throughput: {reads / elapsed:.0f} reads/sec")
    except:
        pass

    return elapsed


if __name__ == "__main__":
    profile_write_remap()
