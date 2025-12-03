#!/usr/bin/env python3
"""
Plot RNA-seq benchmark: WASP1 vs STAR+WASP vs WASP2-Rust
Based on HG00731 sample (56M reads) from Asiimwe & Dobin, bioRxiv 2024
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Published results from Asiimwe & Dobin, bioRxiv 2024
# Sample: HG00731, 56M reads, 8 threads
PUBLISHED_RESULTS = {
    'WASP1': 1541,      # 25:41.13
    'STAR+WASP': 332,   # 5:31.67
    'STAR_only': 296,   # 4:55.91 (reference, no ASE correction)
}

def create_comparison_plot(wasp2_rust_time: float, output_dir: str = 'plots'):
    """Create bar chart comparing methods."""
    Path(output_dir).mkdir(exist_ok=True)

    methods = ['WASP1\n(original)', 'STAR+WASP', 'WASP2-Rust\nv1.3.0']
    times = [PUBLISHED_RESULTS['WASP1'], PUBLISHED_RESULTS['STAR+WASP'], wasp2_rust_time]
    colors = ['gray', 'steelblue', 'crimson']

    fig, ax = plt.subplots(figsize=(10, 7))

    bars = ax.bar(methods, times, color=colors, edgecolor='black', linewidth=1.2)

    # Add time labels on bars
    for bar, time in zip(bars, times):
        height = bar.get_height()
        minutes = int(time // 60)
        seconds = int(time % 60)
        ax.text(bar.get_x() + bar.get_width()/2., height + 30,
                f'{minutes}:{seconds:02d}\n({time:.0f}s)',
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Add speedup annotations
    wasp1_time = PUBLISHED_RESULTS['WASP1']
    for i, (method, time) in enumerate(zip(methods, times)):
        if method != 'WASP1\n(original)':
            speedup = wasp1_time / time
            ax.text(i, time / 2, f'{speedup:.1f}x\nfaster',
                    ha='center', va='center', fontsize=14, fontweight='bold',
                    color='white')

    ax.set_ylabel('Time (seconds)', fontsize=14)
    ax.set_title('RNA-seq ASE Correction Performance\nHG00731 Sample (56M reads, 8 threads)',
                 fontsize=16, fontweight='bold')
    ax.set_ylim(0, max(times) * 1.15)

    # Add horizontal line for STAR-only reference
    ax.axhline(y=PUBLISHED_RESULTS['STAR_only'], color='green', linestyle='--',
               linewidth=2, label=f'STAR-only (no ASE): {PUBLISHED_RESULTS["STAR_only"]}s')
    ax.legend(loc='upper right', fontsize=11)

    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/star_wasp_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/star_wasp_comparison.png")
    plt.close()

    # Print summary
    print("\n" + "=" * 60)
    print("RNA-seq ASE Correction Benchmark Summary")
    print("=" * 60)
    print(f"Sample: HG00731 (56M paired-end reads)")
    print(f"Threads: 8")
    print("-" * 60)
    print(f"{'Method':<20} {'Time':>10} {'vs WASP1':>12}")
    print("-" * 60)
    for method, time in zip(['WASP1', 'STAR+WASP', 'WASP2-Rust'], times):
        speedup = wasp1_time / time
        print(f"{method:<20} {time:>7.0f}s {speedup:>10.1f}x")
    print("=" * 60)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        wasp2_time = float(sys.argv[1])
    else:
        # Placeholder - will update when benchmark completes
        wasp2_time = 200  # Estimate based on scaling data
        print("Using estimated WASP2-Rust time. Pass actual time as argument.")

    create_comparison_plot(wasp2_time)
