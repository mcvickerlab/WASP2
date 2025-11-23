import argparse
import subprocess
import time
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def run_one(mode: str, bam: Path, vcf: Path, bed: Path, out_path: Path, vcf_bed: Optional[Path], intersect_bed: Optional[Path], threads: int = 1) -> tuple[float, Path]:
    """Run counting with given mode ('python' or 'rust') and return wall time + output path."""
    import os
    cmd = [
        "python",
        "-m",
        "src.counting",
        "count-variants",
        str(bam),
        str(vcf),
        "--regions",
        str(bed),
        "--out_file",
        str(out_path),
    ]
    if vcf_bed:
        cmd.extend(["--vcf-bed", str(vcf_bed)])
    if intersect_bed:
        cmd.extend(["--intersect-bed", str(intersect_bed)])
    if mode == "python":
        cmd.append("--no-rust")

    # Set thread count via environment variable
    env = os.environ.copy()
    env["WASP2_RUST_THREADS"] = str(threads)

    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, env=env)
    wall = time.perf_counter() - t0
    return wall, out_path


def main():
    parser = argparse.ArgumentParser(description="Benchmark counting (Python vs Rust) on chr10 test data")
    parser.add_argument("--bam", type=Path, default=ROOT / "test_data" / "CD4_ATACseq_Day1_merged_filtered.sort.bam")
    parser.add_argument("--vcf", type=Path, default=ROOT / "test_data" / "filter_chr10.vcf")
    parser.add_argument("--bed", type=Path, default=ROOT / "test_data" / "NA12878_snps_chr10.bed")
    parser.add_argument("--out-csv", type=Path, default=ROOT / "analysis" / "counting_bench.csv")
    parser.add_argument("--out-png", type=Path, default=ROOT / "analysis" / "counting_bench.png")
    parser.add_argument("--vcf-bed", type=Path, default=None, help="Optional precomputed VCF bed to skip vcf_to_bed")
    parser.add_argument("--intersect-bed", type=Path, default=None, help="Optional precomputed intersect to skip bedtools intersect")
    args = parser.parse_args()

    rows = []
    # Python baseline (threads=1)
    out_file = ROOT / "test_python_counts.tsv"
    wall, _ = run_one("python", args.bam, args.vcf, args.bed, out_file, args.vcf_bed, args.intersect_bed, threads=1)
    rows.append({
        "method": "python_t1",
        "threads": 1,
        "wall_time_sec": wall,
    })

    # Rust with different thread counts
    for threads in [1, 8, 16]:
        out_file = ROOT / f"test_rust_t{threads}_counts.tsv"
        wall, _ = run_one("rust", args.bam, args.vcf, args.bed, out_file, args.vcf_bed, args.intersect_bed, threads=threads)
        rows.append({
            "method": f"rust_t{threads}",
            "threads": threads,
            "wall_time_sec": wall,
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.out_csv, index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["C0", "C1", "C2", "C3"]
    x_labels = ["Python\n(1 thread)", "Rust\n(1 thread)", "Rust\n(8 threads)", "Rust\n(16 threads)"]
    bars = ax.bar(range(len(df)), df["wall_time_sec"], color=colors)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Counting benchmark (chr10 test, with/without threading)")
    for i, v in enumerate(df["wall_time_sec"]):
        ax.text(i, v, f"{v:.2f}s", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_png, dpi=150)
    plt.close(fig)

    print(df.to_string(index=False))
    print(f"\nWrote {args.out_csv} and {args.out_png}")


if __name__ == "__main__":
    main()
