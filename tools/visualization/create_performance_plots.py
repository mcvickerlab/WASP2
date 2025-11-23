#!/usr/bin/env python3
"""
Create Nature-style performance plots for WASP2 Rust optimization
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# Nature journal style settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Colorblind-friendly palette (Wong 2011)
PYTHON_COLOR = '#0173B2'  # Blue
RUST_COLOR = '#DE8F05'    # Orange
COLORS = {
    'python': PYTHON_COLOR,
    'rust': RUST_COLOR,
    'speedup': '#029E73'  # Green
}

# Performance data (from benchmarks)
PERFORMANCE_DATA = {
    'Counting': {
        'python_time': 9.26,      # seconds
        'python_std': 0.5,
        'rust_time': 1.30,
        'rust_std': 0.1,
        'python_mem': 639,        # MB
        'rust_mem': 300,
        'speedup': 7.12
    },
    'Mapping': {
        'python_time': 0.147,     # seconds
        'python_std': 0.015,
        'rust_time': 0.0323,
        'rust_std': 0.005,
        'python_mem': 100,
        'rust_mem': 50,
        'speedup': 4.55
    }
}

def create_performance_plot():
    """Plot 1: Performance comparison (time)"""
    fig, ax = plt.subplots(figsize=(3.5, 2.5), dpi=300)

    stages = list(PERFORMANCE_DATA.keys())
    x = np.arange(len(stages))
    width = 0.35

    # Get data
    python_times = [PERFORMANCE_DATA[s]['python_time'] * 1000 for s in stages]  # Convert to ms
    python_stds = [PERFORMANCE_DATA[s]['python_std'] * 1000 for s in stages]
    rust_times = [PERFORMANCE_DATA[s]['rust_time'] * 1000 for s in stages]
    rust_stds = [PERFORMANCE_DATA[s]['rust_std'] * 1000 for s in stages]

    # Create bars
    bars1 = ax.bar(x - width/2, python_times, width, yerr=python_stds,
                   label='Python', color=COLORS['python'], capsize=3,
                   error_kw={'linewidth': 0.5})
    bars2 = ax.bar(x + width/2, rust_times, width, yerr=rust_stds,
                   label='Rust', color=COLORS['rust'], capsize=3,
                   error_kw={'linewidth': 0.5})

    # Styling
    ax.set_ylabel('Time (ms, log scale)', fontsize=8)
    ax.set_xlabel('Pipeline Stage', fontsize=8)
    ax.set_title('Performance Comparison', fontsize=9, fontweight='bold')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add speedup annotations
    for i, stage in enumerate(stages):
        speedup = PERFORMANCE_DATA[stage]['speedup']
        y_pos = max(python_times[i], rust_times[i]) * 1.5
        ax.text(i, y_pos, f'{speedup:.1f}x', ha='center', va='bottom',
                fontsize=7, fontweight='bold', color=COLORS['speedup'])

    plt.tight_layout()
    return fig

def create_memory_plot():
    """Plot 2: Memory usage comparison"""
    fig, ax = plt.subplots(figsize=(3.5, 2.5), dpi=300)

    stages = list(PERFORMANCE_DATA.keys())
    x = np.arange(len(stages))
    width = 0.35

    # Get data
    python_mem = [PERFORMANCE_DATA[s]['python_mem'] for s in stages]
    rust_mem = [PERFORMANCE_DATA[s]['rust_mem'] for s in stages]

    # Create bars
    bars1 = ax.bar(x - width/2, python_mem, width,
                   label='Python', color=COLORS['python'])
    bars2 = ax.bar(x + width/2, rust_mem, width,
                   label='Rust', color=COLORS['rust'])

    # Styling
    ax.set_ylabel('Memory (MB)', fontsize=8)
    ax.set_xlabel('Pipeline Stage', fontsize=8)
    ax.set_title('Memory Usage Reduction', fontsize=9, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add reduction percentages
    for i, stage in enumerate(stages):
        reduction = (1 - PERFORMANCE_DATA[stage]['rust_mem'] / PERFORMANCE_DATA[stage]['python_mem']) * 100
        y_pos = max(python_mem[i], rust_mem[i]) * 1.05
        ax.text(i, y_pos, f'-{reduction:.0f}%', ha='center', va='bottom',
                fontsize=7, fontweight='bold', color='#CC3311')

    plt.tight_layout()
    return fig

def create_speedup_plot():
    """Plot 3: Speedup fold-change"""
    fig, ax = plt.subplots(figsize=(3.5, 2.5), dpi=300)

    stages = list(PERFORMANCE_DATA.keys())
    speedups = [PERFORMANCE_DATA[s]['speedup'] for s in stages]

    # Create bars
    colors = [COLORS['speedup'] if s > 1 else '#CC3311' for s in speedups]
    bars = ax.bar(stages, speedups, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

    # Reference line at 1x (no improvement)
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='No improvement')

    # Styling
    ax.set_ylabel('Speedup (fold-change)', fontsize=8)
    ax.set_xlabel('Pipeline Stage', fontsize=8)
    ax.set_title('Performance Improvement', fontsize=9, fontweight='bold')
    ax.set_yscale('log', base=2)
    ax.set_ylim(bottom=0.5)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add value labels
    for i, (stage, speedup) in enumerate(zip(stages, speedups)):
        ax.text(i, speedup * 1.1, f'{speedup:.1f}x', ha='center', va='bottom',
                fontsize=8, fontweight='bold')

    plt.tight_layout()
    return fig

def create_scaling_plot():
    """Plot 4: Whole-genome extrapolation"""
    fig, ax = plt.subplots(figsize=(3.5, 2.5), dpi=300)

    # Extrapolate from chr10 test data to whole genome
    # chr10 is ~133 Mbp, human genome is ~3,088 Mbp (23x larger)
    scale_factor = 23

    # Test data point (chr10)
    test_reads = 2409
    python_test_time = (PERFORMANCE_DATA['Counting']['python_time'] +
                        PERFORMANCE_DATA['Mapping']['python_time'])  # ~9.4s
    rust_test_time = (PERFORMANCE_DATA['Counting']['rust_time'] +
                      PERFORMANCE_DATA['Mapping']['rust_time'])  # ~1.3s

    # Extrapolation (linear scaling with genome size)
    read_counts = np.array([test_reads, test_reads * 10, test_reads * 23, test_reads * 50])
    python_times = read_counts / test_reads * python_test_time / 60  # Convert to minutes
    rust_times = read_counts / test_reads * rust_test_time / 60

    # Plot lines
    ax.plot(read_counts, python_times, 'o-', color=COLORS['python'],
            label='Python', linewidth=1.5, markersize=5)
    ax.plot(read_counts, rust_times, 'o-', color=COLORS['rust'],
            label='Rust', linewidth=1.5, markersize=5)

    # Mark actual measurement
    ax.plot(test_reads, python_test_time/60, 'o', color=COLORS['python'],
            markersize=8, markerfacecolor='white', markeredgewidth=2)
    ax.plot(test_reads, rust_test_time/60, 'o', color=COLORS['rust'],
            markersize=8, markerfacecolor='white', markeredgewidth=2)

    # Styling
    ax.set_xlabel('Read Pairs Processed', fontsize=8)
    ax.set_ylabel('Time (minutes)', fontsize=8)
    ax.set_title('Scalability: Whole Genome Processing', fontsize=9, fontweight='bold')
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add annotation for whole genome
    wg_idx = 2  # Index for whole genome point
    ax.annotate('Whole genome\n(~55k pairs)',
                xy=(read_counts[wg_idx], rust_times[wg_idx]),
                xytext=(read_counts[wg_idx] * 0.7, rust_times[wg_idx] * 3),
                fontsize=6, ha='center',
                arrowprops=dict(arrowstyle='->', linewidth=0.5, color='black'))

    # Add speedup annotation
    speedup_text = f'{python_times[wg_idx] / rust_times[wg_idx]:.1f}x faster'
    ax.text(0.95, 0.95, speedup_text, transform=ax.transAxes,
            fontsize=7, fontweight='bold', ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3, linewidth=0.5))

    plt.tight_layout()
    return fig

def main():
    """Generate all plots"""
    output_dir = Path(__file__).parent.parent.parent / 'docs' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Generating Nature-style performance plots...")
    print("=" * 70)

    # Plot 1: Performance
    print("\n1. Performance Comparison (Time)...")
    fig1 = create_performance_plot()
    fig1.savefig(output_dir / 'performance_comparison.pdf', bbox_inches='tight')
    fig1.savefig(output_dir / 'performance_comparison.png', bbox_inches='tight')
    print(f"   ✓ Saved to {output_dir}/performance_comparison.[pdf|png]")

    # Plot 2: Memory
    print("\n2. Memory Usage Reduction...")
    fig2 = create_memory_plot()
    fig2.savefig(output_dir / 'memory_usage.pdf', bbox_inches='tight')
    fig2.savefig(output_dir / 'memory_usage.png', bbox_inches='tight')
    print(f"   ✓ Saved to {output_dir}/memory_usage.[pdf|png]")

    # Plot 3: Speedup
    print("\n3. Performance Improvement (Fold-Change)...")
    fig3 = create_speedup_plot()
    fig3.savefig(output_dir / 'speedup_foldchange.pdf', bbox_inches='tight')
    fig3.savefig(output_dir / 'speedup_foldchange.png', bbox_inches='tight')
    print(f"   ✓ Saved to {output_dir}/speedup_foldchange.[pdf|png]")

    # Plot 4: Scaling
    print("\n4. Scalability (Whole Genome Extrapolation)...")
    fig4 = create_scaling_plot()
    fig4.savefig(output_dir / 'scaling_extrapolation.pdf', bbox_inches='tight')
    fig4.savefig(output_dir / 'scaling_extrapolation.png', bbox_inches='tight')
    print(f"   ✓ Saved to {output_dir}/scaling_extrapolation.[pdf|png]")

    print("\n" + "=" * 70)
    print("✓ All plots generated successfully!")
    print(f"\nOutput directory: {output_dir}")
    print("\nPlot specifications:")
    print("  - Format: PDF + PNG")
    print("  - DPI: 300 (publication quality)")
    print("  - Style: Nature journal standards")
    print("  - Colors: Colorblind-friendly (Wong 2011 palette)")
    print("  - Fonts: Arial, 8-10pt")

if __name__ == '__main__':
    main()
