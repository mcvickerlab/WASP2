#!/usr/bin/env python3
"""Create publication-ready combined benchmark figure"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

ROOT = Path(__file__).parent

# Load benchmark data
counting_df = pd.read_csv(ROOT / "counting_bench.csv")
mapping_df = pd.read_csv(ROOT / "mapping_filter_bench.csv")
analysis_df = pd.read_csv(ROOT / "analysis_bench.csv")

# Create publication figure
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Counting
ax = axes[0]
python_time = counting_df[counting_df['method'] == 'python_t1']['wall_time_sec'].values[0]
rust_t1 = counting_df[counting_df['method'] == 'rust_t1']['wall_time_sec'].values[0]

x_pos = [0, 1]
times = [python_time, rust_t1]
colors = ['#3498db', '#e74c3c']
bars = ax.bar(x_pos, times, color=colors, width=0.6)

ax.set_ylabel('Wall time (s)', fontsize=12, fontweight='bold')
ax.set_title('A. Allele Counting\n(111K SNPs, chr10)', fontsize=13, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(['Python', 'Rust'], fontsize=11)
ax.set_ylim(0, max(times) * 1.15)

# Add speedup annotation
speedup = python_time / rust_t1
ax.text(0.5, max(times) * 0.9, f'{speedup:.1f}× faster',
        ha='center', fontsize=14, weight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.3))

# Add time labels
for i, v in enumerate(times):
    ax.text(i, v, f'{v:.2f}s', ha='center', va='bottom', fontsize=10, weight='bold')

# Panel B: Mapping Filter
ax = axes[1]
python_time = mapping_df[mapping_df['method'] == 'python_t1']['wall_time_sec'].values[0]
rust_t16 = mapping_df[mapping_df['method'] == 'rust_t16']['wall_time_sec'].values[0]

times = [python_time, rust_t16]
bars = ax.bar(x_pos, times, color=colors, width=0.6)

ax.set_ylabel('Wall time (s)', fontsize=12, fontweight='bold')
ax.set_title('B. Mapping Filter\n(20K pairs, 16 threads)', fontsize=13, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(['Python\n(1 thread)', 'Rust\n(16 threads)'], fontsize=11)
ax.set_ylim(0, max(times) * 1.15)

speedup = python_time / rust_t16
ax.text(0.5, max(times) * 0.9, f'{speedup:.1f}× faster',
        ha='center', fontsize=14, weight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.3))

for i, v in enumerate(times):
    ax.text(i, v, f'{v:.2f}s', ha='center', va='bottom', fontsize=10, weight='bold')

# Panel C: Analysis
ax = axes[2]
python_time = analysis_df[analysis_df['method'] == 'python']['wall_time_sec'].values[0]
rust_time = analysis_df[analysis_df['method'] == 'rust']['wall_time_sec'].values[0]

times = [python_time, rust_time]
bars = ax.bar(x_pos, times, color=colors, width=0.6)

ax.set_ylabel('Wall time (s)', fontsize=12, fontweight='bold')
ax.set_title('C. Allelic Imbalance\n(43 regions, β-binomial)', fontsize=13, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(['Python', 'Rust'], fontsize=11)
ax.set_ylim(0, max(times) * 1.15)

speedup = python_time / rust_time
ax.text(0.5, max(times) * 0.9, f'{speedup:.1f}× faster',
        ha='center', fontsize=14, weight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.3))

for i, v in enumerate(times):
    ax.text(i, v, f'{v:.2f}s', ha='center', va='bottom', fontsize=10, weight='bold')

# Overall styling
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

fig.suptitle('WASP2 Rust Acceleration Performance',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save
output_path = ROOT / 'publication_benchmarks.png'
fig.savefig(output_path, dpi=300, bbox_inches='tight')
output_pdf = ROOT / 'publication_benchmarks.pdf'
fig.savefig(output_pdf, bbox_inches='tight')

print(f"Publication figure saved to:")
print(f"  {output_path}")
print(f"  {output_pdf}")

plt.close()
