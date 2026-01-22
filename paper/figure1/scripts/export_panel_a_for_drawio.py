#!/usr/bin/env python3
"""
Export Figure 1 Panel A for DrawIO Import

This script generates standalone exports of the WASP2 architecture diagram
(Panel A from Figure 1) in multiple formats suitable for:
- DrawIO import and manual editing (SVG)
- High-resolution presentation/documentation (PNG at 600 DPI)
- Vector publication format (PDF)

Usage:
    python export_panel_a_for_drawio.py

Output files are saved to:
    paper/figure1/exports/panel_a_standalone.{svg,png,pdf}
"""
import sys
from pathlib import Path

# Add paper directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server use

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Import the panel_a_architecture function and style setup
from figure1.scripts.generate_figure1 import panel_a_architecture, setup_style, C, FONTS


def create_panel_a_standalone():
    """Create a standalone figure with just Panel A (architecture diagram)."""
    setup_style()

    # Figure size optimized for architecture diagram
    # Nature Methods double column width: 180mm = 7.09 inches
    # Panel A aspect ratio is approximately 1:0.55 (width:height)
    fig_width = 7.09  # Double column width
    fig_height = 4.5  # Height that captures the diagram well

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Single axes spanning the figure
    ax = fig.add_subplot(111)

    # Generate the architecture diagram
    panel_a_architecture(ax)

    # Remove the panel label 'a' for standalone export (optional - comment out if you want it)
    # The panel_a_architecture function adds 'a' at position (-0.02, 1.02)
    # We can leave it or remove it by redrawing

    # Adjust layout to minimize whitespace
    plt.tight_layout(pad=0.5)

    return fig


def export_panel_a():
    """Export Panel A in multiple formats for DrawIO and publication use."""
    exports_dir = Path(__file__).parent.parent / 'exports'
    exports_dir.mkdir(parents=True, exist_ok=True)

    print("Generating Panel A standalone figure...")
    fig = create_panel_a_standalone()

    # Export formats and settings
    exports = [
        {
            'path': exports_dir / 'panel_a_standalone.svg',
            'format': 'svg',
            'dpi': 300,  # SVG is vector but DPI affects text rendering
            'description': 'SVG for DrawIO import (vector, editable)'
        },
        {
            'path': exports_dir / 'panel_a_standalone.png',
            'format': 'png',
            'dpi': 600,  # High resolution for maximum quality
            'description': 'High-res PNG (600 DPI raster)'
        },
        {
            'path': exports_dir / 'panel_a_standalone.pdf',
            'format': 'pdf',
            'dpi': 300,  # PDF is vector
            'description': 'PDF vector format'
        },
    ]

    for export in exports:
        print(f"  Exporting: {export['path'].name} - {export['description']}")
        fig.savefig(
            export['path'],
            format=export['format'],
            dpi=export['dpi'],
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            transparent=False,
        )
        print(f"    Saved: {export['path']}")

    plt.close(fig)

    # Print file sizes
    print("\nExport summary:")
    for export in exports:
        size_kb = export['path'].stat().st_size / 1024
        print(f"  {export['path'].name}: {size_kb:.1f} KB")

    return exports_dir


def create_drawio_import_notes():
    """Generate notes for importing into DrawIO."""
    notes = """
DrawIO Import Notes for Panel A
================================

The SVG export is optimized for import into DrawIO (diagrams.net).

Quick Import Steps:
1. Open DrawIO (https://app.diagrams.net or desktop app)
2. File -> Import From -> Device
3. Select panel_a_standalone.svg
4. Click Open

After Import - Editing:
1. Select All (Ctrl+A)
2. Ungroup (Ctrl+Shift+U) - may need to repeat multiple times
3. Individual elements (boxes, text, arrows) are now editable
4. Use Format panel to change colors, fonts, line styles

Tips for Best Results:
- The SVG preserves vector shapes and text
- Colors use the colorblind-safe Paul Tol palette
- Font is Arial (or closest system match)
- Some elements may need ungrouping multiple times

Re-export for Publication:
1. File -> Export As -> PDF
2. Or File -> Export As -> PNG (set scale to 300% for 900 DPI equivalent)

Color Reference (from source):
- Component 1 (Blue): #4477AA / #DCE9F4 (border/fill)
- Component 2 (Green): #228833 / #D9EDDB
- Component 3 (Gold): #CCBB44 / #F5F2DC
- Single-cell (Magenta): #AA3377 / #F4E4ED
- Container: #F5F8FC / #B8C9DC (fill/border)
- Input/Output: #555555 / #E8E8E8
"""
    return notes


if __name__ == '__main__':
    exports_dir = export_panel_a()
    print(f"\nAll exports saved to: {exports_dir}")
    print("\nSee README.md in exports directory for DrawIO import instructions.")
