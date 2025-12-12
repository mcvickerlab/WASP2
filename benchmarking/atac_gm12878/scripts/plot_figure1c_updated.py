#!/usr/bin/env python3
"""
Generate updated Figure 1C: SNV Read Retention by Pipeline
Shows 4 pipelines with their SNV retention rates.
INDEL retention is placeholder for WASP2-Rust (+INDEL) only.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
SCRIPT_DIR = Path(__file__).parent
DATA_FILE = SCRIPT_DIR.parent.parent.parent / "paper/figure1/data/snv_indel_read_counts.tsv"
OUTPUT_DIR = SCRIPT_DIR.parent / "figures"

# Colors matching Panel B
COLORS = {
    'wasp1': '#808080',           # Gray
    'wasp2python': '#4169E1',     # Blue
    'wasp2rust_snp': '#008B8B',   # Teal/Dark Cyan
    'wasp2rust_indel': '#C71585', # Magenta/Pink
}

LABELS = {
    'wasp1': 'WASP1',
    'wasp2python': 'WASP2-Python',
    'wasp2rust_snp': 'WASP2-Rust\n(SNV)',
    'wasp2rust_indel': 'WASP2-Rust\n(+INDEL)',
}

def setup_style():
    """Nature Methods style."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 1,
        'xtick.labelsize': 9,
        'ytick.labelsize': 10,
        'figure.dpi': 150,
        'savefig.dpi': 300,
    })

def load_data():
    """Load SNV retention data from TSV."""
    df = pd.read_csv(DATA_FILE, sep='\t')
    print(f"Loaded data from: {DATA_FILE}")
    print(df)
    return df

def plot_figure1c(df):
    """Generate Figure 1C bar chart."""
    setup_style()

    fig, ax = plt.subplots(figsize=(8, 5))

    # Filter to SNV only
    snv_data = df[df['variant_type'] == 'SNV'].copy()

    # Order pipelines
    pipeline_order = ['wasp1', 'wasp2python', 'wasp2rust_snp', 'wasp2rust_indel']
    snv_data['pipeline'] = pd.Categorical(snv_data['pipeline'], categories=pipeline_order, ordered=True)
    snv_data = snv_data.sort_values('pipeline')

    x = np.arange(len(pipeline_order))
    width = 0.35

    # Plot SNV bars (Original and Retained side by side)
    original_bars = ax.bar(x - width/2, snv_data['reads_pre'] / 1e6, width,
                           label='Original', color='#B0C4DE', edgecolor='black', linewidth=0.5)
    retained_bars = ax.bar(x + width/2, snv_data['reads_post'] / 1e6, width,
                           label='Retained', color='#228B22', edgecolor='black', linewidth=0.5)

    # Add pass rate labels on retained bars
    for i, (bar, rate) in enumerate(zip(retained_bars, snv_data['pass_rate'])):
        height = bar.get_height()
        ax.annotate(f'{rate:.1f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom',
                    fontsize=10, fontweight='bold')

    # Labels and formatting
    ax.set_xlabel('Pipeline', fontsize=12)
    ax.set_ylabel('SNV-overlapping read pairs (millions)', fontsize=12)
    ax.set_title('Figure 1C: SNV Read Retention by Pipeline\n(GM12878 ATAC-seq, 159M reads)', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[p] for p in pipeline_order])
    ax.legend(loc='upper right', frameon=False)

    # Set y-axis limit
    ax.set_ylim(0, 6)

    # Add note about INDEL
    ax.text(0.98, 0.02, 'Note: INDEL retention data pending',
            transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
            style='italic', color='gray')

    plt.tight_layout()

    return fig

def main():
    # Load data
    df = load_data()

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Generate figure
    fig = plot_figure1c(df)

    # Save in multiple formats
    for fmt in ['png', 'pdf', 'svg']:
        output_path = OUTPUT_DIR / f"figure1c_snv_retention.{fmt}"
        fig.savefig(output_path, format=fmt, bbox_inches='tight', dpi=300)
        print(f"Saved: {output_path}")

    plt.close(fig)
    print("\nFigure 1C generation complete!")

if __name__ == "__main__":
    main()
