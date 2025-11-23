import argparse
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
import sys
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from analysis.as_analysis import get_imbalance  # noqa: E402

try:
    import importlib
    import wasp2_rust
    importlib.reload(wasp2_rust)  # Force reload to get latest version
except ImportError:
    wasp2_rust = None


def run_python_analysis(tsv_path: Path, min_count: int, pseudocount: int, method: str):
    """Run Python analysis"""
    t0 = time.perf_counter()
    results = get_imbalance(
        str(tsv_path),
        min_count=min_count,
        pseudocount=pseudocount,
        method=method
    )
    wall = time.perf_counter() - t0
    return results, wall


def run_rust_analysis(tsv_path: Path, min_count: int, pseudocount: int, method: str):
    """Run Rust analysis"""
    if wasp2_rust is None:
        raise RuntimeError("wasp2_rust module not available")

    t0 = time.perf_counter()
    results = wasp2_rust.analyze_imbalance(
        str(tsv_path),
        min_count=min_count,
        pseudocount=pseudocount,
        method=method
    )
    wall = time.perf_counter() - t0
    return results, wall


def main():
    parser = argparse.ArgumentParser(description="Benchmark analysis stage (Python vs Rust)")
    parser.add_argument("--tsv", type=Path, default=ROOT / "test_python_counts.tsv")
    parser.add_argument("--min-count", type=int, default=10, help="Minimum total count threshold")
    parser.add_argument("--pseudocount", type=int, default=1, help="Pseudocount to add")
    parser.add_argument("--method", type=str, default="single", choices=["single", "linear"])
    parser.add_argument("--out-csv", type=Path, default=ROOT / "analysis" / "analysis_bench.csv")
    parser.add_argument("--out-png", type=Path, default=ROOT / "analysis" / "analysis_bench.png")
    args = parser.parse_args()

    if not args.tsv.exists():
        print(f"Error: {args.tsv} not found")
        return

    results_list = []

    # Benchmark Python
    print("Running Python analysis...")
    py_results, py_wall = run_python_analysis(args.tsv, args.min_count, args.pseudocount, args.method)
    print(f"  Python: {py_wall:.2f}s, {len(py_results)} regions analyzed")

    results_list.append({
        "method": "python",
        "wall_time_sec": py_wall,
        "regions_analyzed": len(py_results),
        "min_count": args.min_count,
        "pseudocount": args.pseudocount,
    })

    # Benchmark Rust
    if wasp2_rust is not None:
        print("Running Rust analysis...")
        rust_results, rust_wall = run_rust_analysis(args.tsv, args.min_count, args.pseudocount, args.method)
        print(f"  Rust: {rust_wall:.2f}s, {len(rust_results)} regions analyzed")

        results_list.append({
            "method": "rust",
            "wall_time_sec": rust_wall,
            "regions_analyzed": len(rust_results),
            "min_count": args.min_count,
            "pseudocount": args.pseudocount,
        })
    else:
        print("Skipping Rust (wasp2_rust not available)")

    df = pd.DataFrame(results_list)
    df.to_csv(args.out_csv, index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = ["C0", "C1"]
    x_labels = ["Python", "Rust"] if len(df) == 2 else ["Python"]
    bars = ax.bar(range(len(df)), df["wall_time_sec"], color=colors[:len(df)])
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("Wall time (s)")
    ax.set_title(f"Analysis benchmark ({len(py_results)} regions, method={args.method})")
    for i, wall in enumerate(df["wall_time_sec"]):
        ax.text(i, wall, f"{wall:.2f}s", ha="center", va="bottom", fontsize=10)

    # Add speedup annotation if both methods ran
    if len(df) == 2:
        speedup = df.loc[0, "wall_time_sec"] / df.loc[1, "wall_time_sec"]
        ax.text(0.5, max(df["wall_time_sec"]) * 0.9,
                f"Speedup: {speedup:.1f}×",
                ha="center", fontsize=12, weight="bold",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.3))

    fig.tight_layout()
    fig.savefig(args.out_png, dpi=150)
    plt.close(fig)

    print(f"\n{df.to_string(index=False)}")
    print(f"\nWrote {args.out_csv} and {args.out_png}")


if __name__ == "__main__":
    main()
