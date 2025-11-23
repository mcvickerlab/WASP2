import argparse
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
import sys
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from mapping.filter_remap_reads import filt_remapped_reads  # noqa: E402
from mapping.remap_utils import paired_read_gen  # noqa: E402


def write_bams_from_real(bam_path: Path, tmpdir: Path, n_pairs: int, moved_fraction: float) -> tuple[Path, Path]:
    """Create to_remap/remapped BAMs derived from a real BAM: keep original reads, make remapped copies with WASP names."""
    to_remap_unsorted = tmpdir / "to_remap.unsorted.bam"
    remapped_unsorted = tmpdir / "remapped.unsorted.bam"
    to_remap = tmpdir / "to_remap.bam"
    remapped = tmpdir / "remapped.bam"

    moved_every = max(1, int(1 / moved_fraction)) if moved_fraction > 0 else n_pairs + 1

    with pysam.AlignmentFile(bam_path, "rb") as bam, \
            pysam.AlignmentFile(to_remap_unsorted, "wb", header=bam.header) as out_to, \
            pysam.AlignmentFile(remapped_unsorted, "wb", header=bam.header) as out_remap:
        for i, (r1, r2) in enumerate(paired_read_gen(bam)):
            if i >= n_pairs:
                break

            # write original pair
            out_to.write(r1)
            out_to.write(r2)

            total = 2
            moved = (i % moved_every == 0)
            offset = 5 if moved else 0

            orig_start1 = r1.reference_start
            orig_start2 = r2.reference_start

            for copy_idx in range(1, total + 1):
                r1_copy = pysam.AlignedSegment.fromstring(r1.to_string(), bam.header)
                r2_copy = pysam.AlignedSegment.fromstring(r2.to_string(), bam.header)

                # encode original positions in name
                r1_copy.query_name = f"{r1.query_name}_WASP_{orig_start1}_{orig_start2}_{copy_idx}_{total}"
                r2_copy.query_name = r1_copy.query_name

                start1 = orig_start1 + offset
                start2 = orig_start2 + offset

                r1_copy.reference_start = start1
                r2_copy.reference_start = start2
                r1_copy.next_reference_start = start2
                r2_copy.next_reference_start = start1

                # template lengths
                tlen = max(start1 + r1.query_length, start2 + r2.query_length) - min(start1, start2)
                r1_copy.template_length = tlen
                r2_copy.template_length = -tlen

                out_remap.write(r1_copy)
                out_remap.write(r2_copy)

    # Coordinate-sort then index for pysam.fetch
    pysam.sort("-o", str(to_remap), str(to_remap_unsorted))
    pysam.index(str(to_remap))
    pysam.sort("-o", str(remapped), str(remapped_unsorted))
    pysam.index(str(remapped))

    return to_remap, remapped


def count_pairs(bam_path: Path) -> int:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return bam.count(until_eof=True) // 2


def main():
    parser = argparse.ArgumentParser(description="Benchmark mapping filter (Python vs Rust) on test BAM with synthetic remap copies")
    parser.add_argument("--bam", type=Path, default=ROOT / "test_data" / "CD4_ATACseq_Day1_merged_filtered.sort.bam")
    parser.add_argument("--pairs", type=int, default=20000, help="Number of read pairs to sample")
    parser.add_argument("--moved-fraction", type=float, default=0.2, help="Fraction of pairs forced to move after remap")
    parser.add_argument("--out-csv", type=Path, default=ROOT / "analysis" / "mapping_filter_bench.csv")
    parser.add_argument("--out-png", type=Path, default=ROOT / "analysis" / "mapping_filter_bench.png")
    args = parser.parse_args()

    tmpdir = ROOT / "analysis" / "tmp_mapping_filter"
    tmpdir.mkdir(parents=True, exist_ok=True)
    to_remap, remapped = write_bams_from_real(args.bam, tmpdir, args.pairs, args.moved_fraction)

    results = []

    # Benchmark Python (threads=1)
    out_keep_py = tmpdir / "keep_python_t1.bam"
    t0 = time.perf_counter()
    filt_remapped_reads(str(to_remap), str(remapped), str(out_keep_py), use_rust=False, threads=1)
    wall_py = time.perf_counter() - t0
    kept_py = count_pairs(out_keep_py)
    removed_py = args.pairs - kept_py

    results.append({
        "method": "python_t1",
        "wall_time_sec": wall_py,
        "threads": 1,
        "total_pairs": args.pairs,
        "kept_pairs": kept_py,
        "removed_pairs": removed_py,
        "moved_fraction": args.moved_fraction,
    })

    # Benchmark Rust (threads=1)
    out_keep_rust_t1 = tmpdir / "keep_rust_t1.bam"
    t0 = time.perf_counter()
    filt_remapped_reads(str(to_remap), str(remapped), str(out_keep_rust_t1), use_rust=True, threads=1)
    wall_rust_t1 = time.perf_counter() - t0
    kept_rust_t1 = count_pairs(out_keep_rust_t1)
    removed_rust_t1 = args.pairs - kept_rust_t1

    results.append({
        "method": "rust_t1",
        "wall_time_sec": wall_rust_t1,
        "threads": 1,
        "total_pairs": args.pairs,
        "kept_pairs": kept_rust_t1,
        "removed_pairs": removed_rust_t1,
        "moved_fraction": args.moved_fraction,
    })

    # Benchmark Rust (threads=8)
    out_keep_rust_t8 = tmpdir / "keep_rust_t8.bam"
    t0 = time.perf_counter()
    filt_remapped_reads(str(to_remap), str(remapped), str(out_keep_rust_t8), use_rust=True, threads=8)
    wall_rust_t8 = time.perf_counter() - t0
    kept_rust_t8 = count_pairs(out_keep_rust_t8)
    removed_rust_t8 = args.pairs - kept_rust_t8

    results.append({
        "method": "rust_t8",
        "wall_time_sec": wall_rust_t8,
        "threads": 8,
        "total_pairs": args.pairs,
        "kept_pairs": kept_rust_t8,
        "removed_pairs": removed_rust_t8,
        "moved_fraction": args.moved_fraction,
    })

    # Benchmark Rust (threads=16)
    out_keep_rust_t16 = tmpdir / "keep_rust_t16.bam"
    t0 = time.perf_counter()
    filt_remapped_reads(str(to_remap), str(remapped), str(out_keep_rust_t16), use_rust=True, threads=16)
    wall_rust_t16 = time.perf_counter() - t0
    kept_rust_t16 = count_pairs(out_keep_rust_t16)
    removed_rust_t16 = args.pairs - kept_rust_t16

    results.append({
        "method": "rust_t16",
        "wall_time_sec": wall_rust_t16,
        "threads": 16,
        "total_pairs": args.pairs,
        "kept_pairs": kept_rust_t16,
        "removed_pairs": removed_rust_t16,
        "moved_fraction": args.moved_fraction,
    })

    df = pd.DataFrame(results)
    df.to_csv(args.out_csv, index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["C0", "C1", "C2", "C3"]
    x_labels = ["Python\n(1 thread)", "Rust\n(1 thread)", "Rust\n(8 threads)", "Rust\n(16 threads)"]
    bars = ax.bar(range(len(df)), df["wall_time_sec"], color=colors)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("Wall time (s)")
    ax.set_title(f"Mapping filter benchmark ({args.pairs} pairs, moved {int(args.moved_fraction*100)}%)")
    for i, wall in enumerate(df["wall_time_sec"]):
        ax.text(i, wall, f"{wall:.2f}s", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=150)
    plt.close(fig)

    print(df.to_string(index=False))
    print(f"\nWrote {args.out_csv} and {args.out_png}")


if __name__ == "__main__":
    main()
