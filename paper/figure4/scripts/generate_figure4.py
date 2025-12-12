#!/usr/bin/env python3
"""
Figure 4: Single Cell Analysis

Panel A: Single cell allelic imbalance (placeholder - needs scRNA/scATAC data)
"""
import sys
from pathlib import Path

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from config import (
    COLORS, PLOT_SETTINGS,
    get_plot_path, get_data_path
)


def setup_style():
    """Set up Nature Methods plotting style."""
    plt.rcParams.update({
        'font.family': PLOT_SETTINGS['font_family'],
        'font.sans-serif': PLOT_SETTINGS['font_sans_serif'],
        'font.size': PLOT_SETTINGS['font_size'],
        'axes.labelsize': PLOT_SETTINGS['axes_labelsize'],
        'axes.titlesize': PLOT_SETTINGS['axes_titlesize'],
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': PLOT_SETTINGS['figure_dpi'],
        'savefig.dpi': PLOT_SETTINGS['savefig_dpi'],
    })


def generate_figure4():
    """Generate Figure 4 (placeholder)."""

    setup_style()

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.text(0.5, 0.5, 'Figure 4: Single Cell Analysis\n\n(Data needed: scRNA-seq or scATAC-seq)',
            transform=ax.transAxes, ha='center', va='center',
            fontsize=12, style='italic', color='gray')

    ax.text(0.5, 0.3, 'Potential analyses:\n• Cell-type specific imprinting\n• Pseudobulk aggregation\n• Sparse count handling',
            transform=ax.transAxes, ha='center', va='center',
            fontsize=10, color='gray')

    ax.axis('off')

    fig.suptitle('Figure 4: Single Cell Analysis', fontsize=12, fontweight='bold', y=0.95)

    # Save
    out_dir = get_plot_path(4, 'figure4').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(get_plot_path(4, 'figure4', 'png'), bbox_inches='tight', facecolor='white')
    plt.savefig(get_plot_path(4, 'figure4', 'pdf'), bbox_inches='tight', facecolor='white')

    print(f"Saved: {get_plot_path(4, 'figure4', 'png')}")
    print(f"Saved: {get_plot_path(4, 'figure4', 'pdf')}")

    return fig


if __name__ == '__main__':
    generate_figure4()
