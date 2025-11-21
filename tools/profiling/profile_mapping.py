#!/usr/bin/env python3
"""
Profile WASP2 mapping stage to identify bottlenecks
"""
import cProfile
import pstats
import io
import sys
import time
import tracemalloc
from pathlib import Path
from pstats import SortKey

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.mapping.run_mapping import run_make_remap_reads


def profile_mapping():
    """Profile the mapping stage with detailed timing"""

    # Test data
    bam_file = "test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam"
    vcf_file = "test_data/filter_chr10.vcf"
    sample = "NA12878"
    out_dir = "baselines/mapping"

    print("=" * 80)
    print("WASP2 Mapping Stage Profiling")
    print("=" * 80)
    print(f"BAM: {bam_file}")
    print(f"VCF: {vcf_file}")
    print(f"Sample: {sample}")
    print(f"Output: {out_dir}")
    print()

    # Ensure output directory exists
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Start memory tracking
    tracemalloc.start()
    start_time = time.time()

    # Create profiler
    profiler = cProfile.Profile()

    print("Starting profiling...")
    profiler.enable()

    try:
        # Run the mapping stage
        run_make_remap_reads(
            bam_file=bam_file,
            vcf_file=vcf_file,
            samples=sample,
            is_paired=True,
            is_phased=True,
            out_dir=out_dir,
            out_json=f"{out_dir}/wasp_data.json"
        )
    except Exception as e:
        print(f"Error during mapping: {e}")
        import traceback
        traceback.print_exc()
    finally:
        profiler.disable()

    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Print timing summary
    total_time = end_time - start_time
    print()
    print("=" * 80)
    print("TIMING SUMMARY")
    print("=" * 80)
    print(f"Total time: {total_time:.2f}s")
    print(f"Peak memory: {peak / 1024 / 1024:.1f} MB")
    print()

    # Save detailed profile
    profile_file = "baselines/mapping/profile_stats.txt"
    with open(profile_file, "w") as f:
        # Write summary
        f.write("=" * 80 + "\n")
        f.write("WASP2 Mapping Stage Profile\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total time: {total_time:.2f}s\n")
        f.write(f"Peak memory: {peak / 1024 / 1024:.1f} MB\n")
        f.write("\n")

        # Sort by cumulative time
        f.write("=" * 80 + "\n")
        f.write("TOP FUNCTIONS BY CUMULATIVE TIME\n")
        f.write("=" * 80 + "\n")
        stats = pstats.Stats(profiler, stream=f)
        stats.strip_dirs()
        stats.sort_stats(SortKey.CUMULATIVE)
        stats.print_stats(50)

        f.write("\n" + "=" * 80 + "\n")
        f.write("TOP FUNCTIONS BY TOTAL TIME (self time)\n")
        f.write("=" * 80 + "\n")
        stats.sort_stats(SortKey.TIME)
        stats.print_stats(50)

        f.write("\n" + "=" * 80 + "\n")
        f.write("TOP FUNCTIONS BY CALL COUNT\n")
        f.write("=" * 80 + "\n")
        stats.sort_stats(SortKey.CALLS)
        stats.print_stats(30)

    print(f"Detailed profile saved to: {profile_file}")

    # Print quick summary to console
    print()
    print("=" * 80)
    print("TOP 20 FUNCTIONS (by cumulative time)")
    print("=" * 80)
    s = io.StringIO()
    stats = pstats.Stats(profiler, stream=s)
    stats.strip_dirs()
    stats.sort_stats(SortKey.CUMULATIVE)
    stats.print_stats(20)
    print(s.getvalue())

    return total_time, peak


if __name__ == "__main__":
    profile_mapping()
