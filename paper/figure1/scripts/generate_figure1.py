#!/usr/bin/env python3
"""
Figure 1: WASP2 Software Architecture (v13 - 5-Panel Layout)

Based on comprehensive visual design audit:
- Fixed Component 3 / Output box overlap (was overlapping by 0.07!)
- Reduced component width from 0.22 to 0.20
- Increased inter-component gaps from 0.025 to 0.035
- Fixed single-cell feature box spacing
- Increased sub-step font to 5.5pt for legibility
- Added spacing constants for maintainability
- Proper arrow buffer zones
- FIXED: Logo now positioned INSIDE container box (top-left with proper padding)

v13 Updates:
- Expanded to 5-panel layout (A-E) following STAR+WASP paper convention
- Panel A: Architecture diagram (full width, top)
- Panel B: ATAC-seq runtime benchmark (bottom-left)
- Panel C: ATAC-seq SNV retention (bottom-center-left)
- Panel D: RNA-seq runtime benchmark (bottom-center-right)
- Panel E: RNA-seq SNV retention (bottom-right)

Nature Methods compliant: 180mm width, 300dpi, 5-7pt fonts.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle, ConnectionPatch, Ellipse
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
import pandas as pd
from config import get_plot_path, get_data_path, REPO_ROOT

# Path to WASP2 logo
LOGO_PATH = REPO_ROOT / 'doc/wasp2_hex_logo_v1.png'
FIG1_BENCH_DIR = get_data_path(1, 'benchmarks')


def first_existing(*paths):
    """Return the first existing Path from the given candidates."""
    for p in paths:
        if p is None:
            continue
        p = Path(p)
        if p.exists():
            return p
    return None

# =========================================
# SPACING CONSTANTS (audit recommendations)
# =========================================
SPACING = {
    'component_gap': 0.035,      # Between major components (was 0.025)
    'element_padding': 0.02,     # Inside boxes to text
    'arrow_buffer': 0.015,       # Gap before arrow touches element (was 0.005)
    'section_gap': 0.04,         # Between pipeline and SC strip
    'text_margin': 0.015,        # Text to any edge
}

FONTS = {
    'logo_subtitle': 7,
    'component_title': 7,
    'sub_step': 5.5,             # Increased from 5 for legibility
    'feature_text': 5,
    'sc_title': 6,
    'sc_feature_title': 5,
    'sc_feature_desc': 4,
}


def setup_style():
    """Nature Methods style."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 6,
        'axes.labelsize': 6,
        'axes.titlesize': 7,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.5,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'legend.fontsize': 5,
        'figure.dpi': 150,
        'savefig.dpi': 300,
    })


# Colorblind-safe palette (Paul Tol bright)
C = {
    # Main components
    'comp1': '#4477AA',         # Blue - Read Filtering
    'comp1_bg': '#DCE9F4',
    'comp2': '#228833',         # Green - Counting
    'comp2_bg': '#D9EDDB',
    'comp3': '#CCBB44',         # Yellow/Gold - AI Detection
    'comp3_bg': '#F5F2DC',

    # Single-cell
    'singlecell': '#AA3377',    # Magenta
    'singlecell_bg': '#F4E4ED',

    # Container
    'container_bg': '#F5F8FC',
    'container_border': '#B8C9DC',

    # Input/Output
    'input': '#555555',
    'input_bg': '#E8E8E8',
    'output': '#555555',

    # UI elements
    'text': '#222222',
    'text_light': '#555555',
    'arrow': '#666666',
    'white': '#FFFFFF',

    # Benchmark
    'wasp1': '#888888',
    'rust_snv': '#E69F00',       # Orange - distinct from Python blue
}


def draw_rounded_box(ax, x, y, w, h, facecolor, edgecolor=None, lw=1.2, alpha=1.0, zorder=1):
    """Draw a rounded rectangle."""
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        facecolor=facecolor,
        edgecolor=edgecolor if edgecolor else facecolor,
        linewidth=lw,
        alpha=alpha,
        zorder=zorder
    )
    ax.add_patch(box)
    return box


def draw_number_badge(ax, x, y, num, color, radius=0.022):
    """Draw a numbered circle badge."""
    circle = Circle((x, y), radius, facecolor=color, edgecolor='white', lw=1.5, zorder=10)
    ax.add_patch(circle)
    ax.text(x, y, str(num), ha='center', va='center',
            fontsize=7, fontweight='bold', color='white', zorder=11)


def draw_arrow(ax, start, end, color='#666666', lw=1.2, style='->', shrinkA=0, shrinkB=0):
    """Draw an arrow with optional shrink for buffer zones."""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                               shrinkA=shrinkA, shrinkB=shrinkB, mutation_scale=12))


def panel_a_architecture(ax):
    """Panel A: Unified WASP2 architecture with pipeline as main flow."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Small vertical nudge for overall layout tightness
    y_shift = 0.015
    # Extra lift for the main pipeline to create more whitespace between
    # stage arrows and the per-stage outputs
    pipeline_shift = 0.02

    # =========================================
    # WASP2 LOGO - Positioned INSIDE container (top-left with padding)
    # =========================================
    logo_x = 0.08   # Shifted right
    logo_y = 0.85 + y_shift + pipeline_shift   # Nudge up slightly
    if LOGO_PATH.exists():
        logo_img = mpimg.imread(str(LOGO_PATH))
        imagebox = OffsetImage(logo_img, zoom=0.055)  # Slightly smaller for header
        ab = AnnotationBbox(imagebox, (logo_x, logo_y),
                           frameon=False, box_alignment=(0.5, 0.5), zorder=5)
        ax.add_artist(ab)
        # Subtitle aligned with logo
        ax.text(logo_x + 0.05, logo_y, 'WASP2', ha='left', va='center',
                fontsize=9, fontweight='bold', color=C['comp1'])
        ax.text(logo_x + 0.135, logo_y, 'Allele-Specific Analysis Toolkit',
                ha='left', va='center', fontsize=FONTS['logo_subtitle'], color=C['text_light'])
    else:
        ax.text(logo_x, logo_y, 'WASP2', ha='left', va='center',
                fontsize=12, fontweight='bold', color=C['comp1'])
        ax.text(logo_x + 0.105, logo_y, 'Allele-Specific Analysis Toolkit',
                ha='left', va='center', fontsize=FONTS['logo_subtitle'], color=C['text_light'])

    # =========================================
    # MAIN CONTAINER BOX (logo inside at top-left)
    # =========================================
    container_x = 0.02
    container_y = 0.08  # Raised up to reduce bottom whitespace
    container_w = 0.96
    container_h = 0.91  # Top at y=0.99, reduced height

    draw_rounded_box(ax, container_x, container_y, container_w, container_h,
                     C['container_bg'], C['container_border'], 1.5, alpha=0.6, zorder=0)

    # =========================================
    # INPUT SECTION (left side, vertical stack with spacing)
    # =========================================
    input_x = 0.04  # Shifted left
    input_w = 0.09
    input_h = 0.06
    input_center_y = 0.50 + y_shift + pipeline_shift  # Nudge up with main pipeline

    # BAM box
    bam_y = input_center_y + 0.08
    draw_rounded_box(ax, input_x, bam_y, input_w, input_h, C['input_bg'], C['input'], 1.0)
    ax.text(input_x + input_w/2, bam_y + input_h/2, 'BAM', ha='center', va='center',
            fontsize=6, fontweight='bold', color=C['text'])

    # VCF box (more space below BAM)
    vcf_y = input_center_y - 0.08
    draw_rounded_box(ax, input_x, vcf_y, input_w, input_h, C['input_bg'], C['input'], 1.0)
    ax.text(input_x + input_w/2, vcf_y + input_h/2, 'VCF/BCF', ha='center', va='center',
            fontsize=6, fontweight='bold', color=C['text'])

    # =========================================
    # MAIN PIPELINE: 3 COMPONENT BOXES (horizontal flow)
    # Fixed spacing to prevent overlap with outputs
    # =========================================
    comp_y = 0.20 + y_shift + pipeline_shift  # Nudge up for overall balance
    comp_h = 0.58  # Taller boxes for more sub-step space
    comp_w = 0.19                           # Slightly narrower for arrow space
    gap = 0.05                              # Larger gap for arrows

    comp1_x = 0.18                           # Shifted left
    comp2_x = comp1_x + comp_w + gap         # = 0.44
    comp3_x = comp2_x + comp_w + gap         # = 0.68

    # Horizontal arrows from BAM and VCF to Component 1 (exact same size as 1→2 and 2→3)
    arrow_start_x = input_x + input_w + 0.015  # Same gap as other arrows
    arrow_end_x = comp1_x - 0.015              # Same gap as other arrows

    # Arrow from BAM - horizontal from center of BAM box
    bam_arrow_y = bam_y + input_h/2
    draw_arrow(ax, (arrow_start_x, bam_arrow_y), (arrow_end_x, bam_arrow_y), C['arrow'], 2.0)
    # Arrow from VCF - horizontal from center of VCF box
    vcf_arrow_y = vcf_y + input_h/2
    draw_arrow(ax, (arrow_start_x, vcf_arrow_y), (arrow_end_x, vcf_arrow_y), C['arrow'], 2.0)

    # --- Component 1: Read Filtering ---
    draw_rounded_box(ax, comp1_x, comp_y, comp_w, comp_h, C['comp1_bg'], C['comp1'], 1.5)
    draw_number_badge(ax, comp1_x + 0.025, comp_y + comp_h - 0.035, 1, C['comp1'])
    ax.text(comp1_x + 0.055, comp_y + comp_h - 0.035, 'Read Filtering', ha='left', va='center',
            fontsize=FONTS['component_title'], fontweight='bold', color=C['comp1'])

    # Internal sub-steps (mini vertical pipeline)
    substep_x = comp1_x + 0.015
    substep_w = comp_w - 0.03
    substep_h = 0.055
    substeps1 = ['Identify Variant Reads', 'Swap Alleles', 'Remap Sequences', 'Concordance Filter']
    substep_start_y = comp_y + comp_h - 0.14  # More space below title

    for i, step in enumerate(substeps1):
        sy = substep_start_y - i * (substep_h + 0.025)  # Match spacing of other boxes
        draw_rounded_box(ax, substep_x, sy, substep_w, substep_h, C['comp1'], None, 0, 0.8)
        ax.text(substep_x + substep_w/2, sy + substep_h/2, step, ha='center', va='center',
                fontsize=FONTS['sub_step'], color=C['white'], fontweight='bold')
        # Mini arrow between substeps (centered in the gap)
        if i < len(substeps1) - 1:
            gap_h = 0.025
            arrow_top = sy - 0.004
            arrow_bottom = sy - gap_h + 0.004
            ax.annotate(
                '',
                xy=(substep_x + substep_w/2, arrow_bottom),
                xytext=(substep_x + substep_w/2, arrow_top),
                arrowprops=dict(arrowstyle='->', color=C['comp1'], lw=0.8, mutation_scale=8),
            )

    # Key feature
    ax.text(comp1_x + comp_w/2, comp_y + 0.025, 'Any aligner supported',
            ha='center', va='center', fontsize=FONTS['feature_text'], fontstyle='italic', color=C['text_light'])

    # Arrow 1→2 - aligned with middle sub-steps, with gaps on both ends
    arrow_y = comp_y + comp_h * 0.55  # Slightly above center, aligned with sub-steps
    draw_arrow(ax, (comp1_x + comp_w + 0.015, arrow_y),
               (comp2_x - 0.015, arrow_y), C['arrow'], 2.0)

    # --- Component 2: Variant Counting ---
    draw_rounded_box(ax, comp2_x, comp_y, comp_w, comp_h, C['comp2_bg'], C['comp2'], 1.5)
    draw_number_badge(ax, comp2_x + 0.025, comp_y + comp_h - 0.035, 2, C['comp2'])
    ax.text(comp2_x + 0.055, comp_y + comp_h - 0.035, 'Variant Counting', ha='left', va='center',
            fontsize=FONTS['component_title'], fontweight='bold', color=C['comp2'])

    # Internal sub-steps
    substeps2 = ['Filter Het Variants', 'Feature Intersection', 'Allele Counting']
    substep_start_y = comp_y + comp_h - 0.14  # More space below title

    for i, step in enumerate(substeps2):
        sy = substep_start_y - i * (substep_h + 0.025)
        draw_rounded_box(ax, comp2_x + 0.015, sy, substep_w, substep_h, C['comp2'], None, 0, 0.8)
        ax.text(comp2_x + 0.015 + substep_w/2, sy + substep_h/2, step, ha='center', va='center',
                fontsize=FONTS['sub_step'], color=C['white'], fontweight='bold')
        if i < len(substeps2) - 1:
            gap_h = 0.025
            arrow_top = sy - 0.004
            arrow_bottom = sy - gap_h + 0.004
            ax.annotate(
                '',
                xy=(comp2_x + 0.015 + substep_w/2, arrow_bottom),
                xytext=(comp2_x + 0.015 + substep_w/2, arrow_top),
                arrowprops=dict(arrowstyle='->', color=C['comp2'], lw=0.8, mutation_scale=8),
            )

    # Key features
    ax.text(comp2_x + comp_w/2, comp_y + 0.05, 'Per-SNV & feature-level',
            ha='center', va='center', fontsize=FONTS['feature_text'], fontstyle='italic', color=C['text_light'])
    ax.text(comp2_x + comp_w/2, comp_y + 0.025, 'Phased haplotype support',
            ha='center', va='center', fontsize=FONTS['feature_text'], fontstyle='italic', color=C['text_light'])

    # Arrow 2→3 - aligned with middle sub-steps, with gaps on both ends
    draw_arrow(ax, (comp2_x + comp_w + 0.015, arrow_y),
               (comp3_x - 0.015, arrow_y), C['arrow'], 2.0)

    # --- Component 3: AI Detection ---
    draw_rounded_box(ax, comp3_x, comp_y, comp_w, comp_h, C['comp3_bg'], C['comp3'], 1.5)
    draw_number_badge(ax, comp3_x + 0.025, comp_y + comp_h - 0.035, 3, C['comp3'])
    ax.text(comp3_x + 0.055, comp_y + comp_h - 0.035, 'AI Detection', ha='left', va='center',
            fontsize=FONTS['component_title'], fontweight='bold', color=C['comp3'])

    # Internal sub-steps
    substeps3 = ['Aggregate Counts', 'Beta-Binomial Model', 'FDR Correction']
    substep_start_y = comp_y + comp_h - 0.14  # More space below title

    for i, step in enumerate(substeps3):
        sy = substep_start_y - i * (substep_h + 0.025)
        draw_rounded_box(ax, comp3_x + 0.015, sy, substep_w, substep_h, C['comp3'], None, 0, 0.8)
        ax.text(comp3_x + 0.015 + substep_w/2, sy + substep_h/2, step, ha='center', va='center',
                fontsize=FONTS['sub_step'], color=C['white'], fontweight='bold')
        if i < len(substeps3) - 1:
            gap_h = 0.025
            arrow_top = sy - 0.004
            arrow_bottom = sy - gap_h + 0.004
            ax.annotate(
                '',
                xy=(comp3_x + 0.015 + substep_w/2, arrow_bottom),
                xytext=(comp3_x + 0.015 + substep_w/2, arrow_top),
                arrowprops=dict(arrowstyle='->', color=C['comp3'], lw=0.8, mutation_scale=8),
            )

    # Key features
    ax.text(comp3_x + comp_w/2, comp_y + 0.05, 'Multi-SNV aggregation',
            ha='center', va='center', fontsize=FONTS['feature_text'], fontstyle='italic', color=C['text_light'])
    ax.text(comp3_x + comp_w/2, comp_y + 0.025, 'Likelihood ratio test',
            ha='center', va='center', fontsize=FONTS['feature_text'], fontstyle='italic', color=C['text_light'])

    # =========================================
    # OUTPUTS (below each component)
    # PI feedback: show outputs per-stage, not stacked at far right.
    # =========================================
    output_w = 0.075
    output_h = 0.055
    # Place outputs above the container bottom, shifted up with the pipeline
    # so stage arrows keep the same length.
    out_y = container_y + 0.01 + y_shift + pipeline_shift

    # Centers aligned to components
    out1_x = comp1_x + comp_w / 2 - output_w / 2
    out2_x = comp2_x + comp_w / 2 - output_w / 2
    out3_x = comp3_x + comp_w / 2 - output_w / 2

    # Downward arrows from each component to its output
    # Shifted down: start further from comp boxes, end closer to output boxes
    # Line width 2.0 to match horizontal arrows
    draw_arrow(
        ax,
        (comp1_x + comp_w / 2, comp_y - 0.035),
        (comp1_x + comp_w / 2, out_y + output_h + 0.010),
        C['arrow'],
        2.0,
    )
    draw_arrow(
        ax,
        (comp2_x + comp_w / 2, comp_y - 0.035),
        (comp2_x + comp_w / 2, out_y + output_h + 0.010),
        C['arrow'],
        2.0,
    )
    draw_arrow(
        ax,
        (comp3_x + comp_w / 2, comp_y - 0.035),
        (comp3_x + comp_w / 2, out_y + output_h + 0.010),
        C['arrow'],
        2.0,
    )

    # Filtered BAM output box (BLUE - from Component 1)
    draw_rounded_box(ax, out1_x, out_y, output_w, output_h, C['comp1_bg'], C['comp1'], 1.0)
    ax.text(out1_x + output_w/2, out_y + output_h/2, 'Filtered\nBAM', ha='center', va='center',
            fontsize=5, fontweight='bold', color=C['text'], linespacing=0.85)

    # Allele-specific Counts output box (GREEN - from Component 2)
    draw_rounded_box(ax, out2_x, out_y, output_w, output_h, C['comp2_bg'], C['comp2'], 1.0)
    ax.text(out2_x + output_w/2, out_y + output_h/2, 'AS\nCounts', ha='center', va='center',
            fontsize=5, fontweight='bold', color=C['text'], linespacing=0.85)

    # AI Results output box (YELLOW - from Component 3)
    draw_rounded_box(ax, out3_x, out_y, output_w, output_h, C['comp3_bg'], C['comp3'], 1.0)
    ax.text(out3_x + output_w/2, out_y + output_h/2, 'AI\nResults', ha='center', va='center',
            fontsize=5, fontweight='bold', color=C['text'], linespacing=0.85)

    # Panel label
    ax.text(-0.02, 1.02, 'a', transform=ax.transAxes, fontsize=8, fontweight='bold')


def panel_b_benchmark(ax):
    """Panel B: Runtime benchmark with actual data from unified pipeline."""
    # Publication source data (tracked)
    wasp1_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp1_perf_log.tsv',
        REPO_ROOT / 'benchmarking/data/aho_perf_logs/wasp1_perf_logs_snp2h5_7881649-Copy1.txt',
    )
    wasp2_py_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp2python_perf_log.tsv',
        REPO_ROOT / 'benchmarking/data/aho_perf_logs/wasp2_perf_logs_7908428-Copy1.txt',
    )
    rust_snp_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp2rust_snv_scaling.tsv',
        REPO_ROOT / 'benchmarking/results/wasp2_unified_scaling_COMBINED.tsv',
    )
    rust_indel_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp2rust_indel_scaling.tsv',
        REPO_ROOT / 'benchmarking/results/wasp2_unified_indel_scaling_COMBINED.tsv',
    )

    if all([wasp1_file, wasp2_py_file, rust_snp_file, rust_indel_file]):
        # WASP1
        wasp1 = pd.read_csv(wasp1_file, sep='\t',
                           names=['n_reads', 'seed', 'total', 'snp2h5', 'intersect', 'remap', 'filter', 'extra'])
        wasp1_g = wasp1.groupby('n_reads')['total'].agg(['mean', 'std']).reset_index()
        wasp1_g['n_reads_m'] = wasp1_g['n_reads'] / 1e6
        wasp1_g['mean_smooth'] = wasp1_g['mean'].rolling(3, center=True, min_periods=1).mean()

        # WASP2-Python
        wasp2_py = pd.read_csv(wasp2_py_file, sep='\t',
                              names=['n_reads', 'seed', 'total', 'intersect', 'remap', 'filter'])
        wasp2_py_g = wasp2_py.groupby('n_reads')['total'].agg(['mean', 'std']).reset_index()
        wasp2_py_g['n_reads_m'] = wasp2_py_g['n_reads'] / 1e6

        # WASP2-Rust SNV (NEW unified pipeline)
        rust_snp = pd.read_csv(rust_snp_file, sep='\t')
        rust_snp_g = rust_snp.groupby('n_reads')['total_s'].agg(['mean', 'std']).reset_index()
        rust_snp_g['n_reads_m'] = rust_snp_g['n_reads'] / 1e6

        # WASP2-Rust SNV+INDEL (NEW unified pipeline with same_locus_slop=10)
        rust_indel = pd.read_csv(rust_indel_file, sep='\t')
        rust_indel_g = rust_indel.groupby('n_reads')['total_s'].agg(['mean', 'std']).reset_index()
        rust_indel_g['n_reads_m'] = rust_indel_g['n_reads'] / 1e6
        rust_indel_g['mean_smooth'] = rust_indel_g['mean'].rolling(3, center=True, min_periods=1).mean()

        # Plot WASP1 (gray, slowest)
        ax.plot(wasp1_g['n_reads_m'], wasp1_g['mean_smooth'], 'o-', color=C['wasp1'],
                label='WASP1', lw=1.5, ms=2.5)
        ax.fill_between(wasp1_g['n_reads_m'], wasp1_g['mean']-wasp1_g['std'],
                       wasp1_g['mean']+wasp1_g['std'], color=C['wasp1'], alpha=0.15)

        # Plot WASP2-Python (Python blue)
        python_blue = '#3776AB'
        ax.plot(wasp2_py_g['n_reads_m'], wasp2_py_g['mean'], 'd-', color=python_blue,
                label='WASP2-Python', lw=1.5, ms=2.5)
        ax.fill_between(wasp2_py_g['n_reads_m'], wasp2_py_g['mean']-wasp2_py_g['std'],
                       wasp2_py_g['mean']+wasp2_py_g['std'], color=python_blue, alpha=0.15)

        # Plot WASP2-Rust SNV (orange - distinct from Python blue)
        ax.plot(rust_snp_g['n_reads_m'], rust_snp_g['mean'], '^-', color=C['rust_snv'],
                label='WASP2-Rust (SNV)', lw=1.5, ms=2.5)

        # Plot WASP2-Rust SNV+INDEL (magenta)
        ax.plot(rust_indel_g['n_reads_m'], rust_indel_g['mean_smooth'], 's-', color=C['singlecell'],
                label='WASP2-Rust (+INDEL)', lw=1.5, ms=2.5)

        ax.set_xlabel('Reads (millions)', fontsize=6)
        ax.set_ylabel('Time (s)', fontsize=6)
        ax.legend(loc='upper left', fontsize=4.5, framealpha=0.95)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.15, lw=0.3)
        ax.set_title('ATAC-seq (GM12878)', fontsize=6, fontweight='bold')

    ax.text(-0.15, 1.05, 'b', transform=ax.transAxes, fontsize=8, fontweight='bold')


def panel_c_filtering(ax):
    """Panel C: ATAC-seq retention of variant-overlapping reads.

    Definition (locked for comparability):
    - "Original" = baseline-filtered variant-overlapping reads that enter the WASP remap+filter
      path (after proper-pair/primary/mapped filters). For INDEL runs this comes directly from
      unified overlap stats; for SNP-only runs we backfill the same denominator using unified stats.
    - "Retained" = reads from that same population that pass WASP remap filtering (remap_keep).
    Units: reads (both mates), stacked by overlap type of the ORIGINAL fragment:
      SNV-only / INDEL-only / SNV+INDEL.
    """
    import json
    from matplotlib.patches import Patch

    type_colors = {
        'snv_only': '#56B4E9',   # Okabe-Ito sky blue
        'indel_only': '#CC79A7', # Okabe-Ito reddish purple
        'both': '#009E73',       # Okabe-Ito bluish green
    }

    def load_split_retention(json_path):
        """Return ({pre_by_type}, {post_by_type}) in READS, or None."""
        if not json_path.exists():
            return None
        d = json.loads(json_path.read_text())

        # New split fields (reads)
        if 'snv_only_overlap_reads_pre' in d:
            pre = {
                'snv_only': float(d.get('snv_only_overlap_reads_pre', 0)),
                'indel_only': float(d.get('indel_only_overlap_reads_pre', 0)),
                'both': float(d.get('snv_indel_overlap_reads_pre', 0)),
            }
            post = {
                'snv_only': float(d.get('snv_only_overlap_reads_post', 0)),
                'indel_only': float(d.get('indel_only_overlap_reads_post', 0)),
                'both': float(d.get('snv_indel_overlap_reads_post', 0)),
            }

            # OLD-STYLE denominator alignment:
            # Use raw overlap population (original_reads - keep_reads) as "Original" total.
            # If split pre total differs (due to baseline filters), rescale pre stacks so that:
            #   sum(pre_by_type) == original_reads - keep_reads.
            if 'original_reads' in d and 'keep_reads' in d:
                old_pre_total = float(d['original_reads'] - d['keep_reads'])
                split_pre_total = sum(pre.values())
                if old_pre_total > 0 and split_pre_total > 0:
                    frac_diff = abs(old_pre_total - split_pre_total) / old_pre_total
                    if frac_diff > 0.005:  # >0.5% mismatch -> rescale
                        scale = old_pre_total / split_pre_total
                        pre = {k: v * scale for k, v in pre.items()}

            return pre, post

        # Fallback for SNP-only JSONs
        if 'original_reads' in d and 'keep_reads' in d and 'remap_keep_reads' in d:
            pre_total = float(d['original_reads'] - d['keep_reads'])
            post_total = float(d['remap_keep_reads'])
            pre = {'snv_only': pre_total, 'indel_only': 0.0, 'both': 0.0}
            post = {'snv_only': post_total, 'indel_only': 0.0, 'both': 0.0}
            return pre, post

        return None

    results = {}

    # WASP2-Python ATAC (fixed SNP-only run)
    py_file = first_existing(
        FIG1_BENCH_DIR / 'panel_c_atac_wasp2python_snv.json',
        REPO_ROOT
        / 'benchmarking/atac_gm12878/results/wasp2python_snp_fixed_2025-12-07_13-12-00/benchmark_results.json',
    )
    r = load_split_retention(py_file) if py_file else None
    if r:
        results['wasp2python'] = {'pre': r[0], 'post': r[1]}

    # WASP2-Rust SNP ATAC (fixed SNP-only run)
    rust_snp_file = first_existing(
        FIG1_BENCH_DIR / 'panel_c_atac_wasp2rust_snv.json',
        REPO_ROOT
        / 'benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-10_00-45-32/benchmark_results.json',
    )
    r = load_split_retention(rust_snp_file) if rust_snp_file else None
    if r:
        results['wasp2rust_snp'] = {'pre': r[0], 'post': r[1]}

    # WASP2-Rust INDEL ATAC (latest fixed run with split fields)
    indel_file = first_existing(FIG1_BENCH_DIR / 'panel_c_atac_wasp2rust_indel.json')
    if indel_file is None:
        rust_indel_files = sorted(
            (REPO_ROOT / 'benchmarking/atac_gm12878/results').glob('wasp2rust_indel_fixed_*/benchmark_results.json')
        )
        indel_file = rust_indel_files[-1] if rust_indel_files else None
    if indel_file is not None:
        r = load_split_retention(indel_file)
        if r:
            results['wasp2rust_indel'] = {'pre': r[0], 'post': r[1]}

    if len(results) >= 2:
        pipeline_order = ['wasp2python', 'wasp2rust_snp']
        pipeline_labels = ['WASP2-\nPython', 'WASP2-Rust\n(SNV)']
        if 'wasp2rust_indel' in results:
            pipeline_order.append('wasp2rust_indel')
            pipeline_labels.append('WASP2-Rust\n(+INDEL)')

        x = np.arange(len(pipeline_order))
        width = 0.35

        pre_totals = np.array(
            [sum(results[p]['pre'].values()) for p in pipeline_order], dtype=float
        )
        post_totals = np.array(
            [sum(results[p]['post'].values()) for p in pipeline_order], dtype=float
        )
        pass_rates = np.where(pre_totals > 0, 100.0 * post_totals / pre_totals, 0.0)

        bottom_pre = np.zeros(len(pipeline_order))
        bottom_post = np.zeros(len(pipeline_order))
        for key in ['snv_only', 'indel_only', 'both']:
            vals_pre = np.array([results[p]['pre'][key] / 1e6 for p in pipeline_order])
            vals_post = np.array([results[p]['post'][key] / 1e6 for p in pipeline_order])
            # Original (pre-filter): light stacks (no hatch per segment)
            ax.bar(
                x - width / 2,
                vals_pre,
                width,
                bottom=bottom_pre,
                color=type_colors[key],
                alpha=0.35,
                edgecolor='none',
            )
            # Retained (post-filter): solid stacks
            ax.bar(
                x + width / 2,
                vals_post,
                width,
                bottom=bottom_post,
                color=type_colors[key],
                alpha=0.95,
                edgecolor='none',
            )
            bottom_pre += vals_pre
            bottom_post += vals_post

        # Overlay uniform hatch across the full Original bar (cleaner than per-segment hatching)
        ax.bar(
            x - width / 2,
            pre_totals / 1e6,
            width,
            bottom=0,
            facecolor='none',
            edgecolor=C['text_light'],
            hatch='///',
            linewidth=0.0,
            zorder=3,
        )

        for i, (post_val, rate) in enumerate(zip(post_totals / 1e6, pass_rates)):
            ax.text(i + width/2, post_val + 0.15, f"{rate:.1f}%",
                    ha='center', fontsize=4.5, color=C['text'], fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels(pipeline_labels, fontsize=4.5)
        ax.set_ylabel('Variant-overlapping\nreads (M)', fontsize=5)
        type_handles = [
            Patch(facecolor=type_colors['snv_only'], label='SNV-only'),
            Patch(facecolor=type_colors['indel_only'], label='INDEL-only'),
            Patch(facecolor=type_colors['both'], label='SNV+INDEL'),
        ]
        ax.legend(handles=type_handles, fontsize=3.8, loc='upper right', framealpha=0.95)
        ax.text(
            0.02,
            0.98,
            'Hatched = Original\nSolid = Retained',
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=3.6,
            color=C['text_light'],
        )
        ax.set_ylim(bottom=0, top=max(pre_totals / 1e6) * 1.25)
        ax.set_title('ATAC-seq (GM12878)', fontsize=6, fontweight='bold')

    ax.text(-0.15, 1.05, 'c', transform=ax.transAxes, fontsize=8, fontweight='bold')


def panel_d_rnaseq_benchmark(ax):
    """Panel D: RNA-seq runtime benchmark - reads from actual benchmark JSONs."""
    import json

    # Load actual benchmark results from JSON files
    times = {}

    # STAR+WASP
    star_wasp_file = first_existing(
        FIG1_BENCH_DIR / 'panel_de_rnaseq_star_wasp.json',
        REPO_ROOT / 'benchmarking/star_wasp_comparison/results/star_wasp_2025-12-04_00-58-54/benchmark_results.json',
    )
    if star_wasp_file is not None and star_wasp_file.exists():
        d = json.loads(star_wasp_file.read_text())
        times['star_wasp'] = d['wall_clock_s']

    # WASP2-Python (from dev branch with multithread support)
    py_file = first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2python.json')
    if py_file is None:
        py_fixed_files = list(
            (REPO_ROOT / 'benchmarking/star_wasp_comparison/results').glob('wasp2python_rnaseq_fixed_*/benchmark_results.json')
        )
        py_file = py_fixed_files[-1] if py_fixed_files else None
    if py_file is not None and py_file.exists():
        d = json.loads(py_file.read_text())
        times['wasp2python'] = d['total_s']

    # WASP2-Rust SNP FAIR (with LoadAndKeep genome caching)
    # Use specific SNP-only run (not the broken INDEL run which lacks indel_mode=True)
    rust_snp_fair_file = first_existing(
        FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_snv.json',
        REPO_ROOT / 'benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-10_19-23-02/benchmark_results.json',
    )
    if rust_snp_fair_file is not None and rust_snp_fair_file.exists():
        d = json.loads(rust_snp_fair_file.read_text())
        times['wasp2rust_snv'] = d['total_s']

    # WASP2-Rust INDEL RNA-seq (v2). Include if old-style retention looks sane (avoid old broken runs).
    indel_file = first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_indel.json')
    if indel_file is None:
        rust_indel_files = list(
            (REPO_ROOT / 'benchmarking/star_wasp_comparison/results').glob('wasp2rust_indel_rnaseq_v2_*/benchmark_results.json')
        )
        indel_file = sorted(rust_indel_files)[-1] if rust_indel_files else None
    if indel_file is not None and indel_file.exists():
        d = json.loads(indel_file.read_text())
        indel_overall = d.get('overall_pass_rate_percent') or d.get('pass_rate_percent')
        if 'snv_only_overlap_pairs_pre' in d and 'snv_only_overlap_pairs_post' in d:
            pre_total = (
                float(d.get('snv_only_overlap_pairs_pre', 0))
                + float(d.get('indel_only_overlap_pairs_pre', 0))
                + float(d.get('snv_indel_overlap_pairs_pre', 0))
            )
            post_total = (
                float(d.get('snv_only_overlap_pairs_post', 0))
                + float(d.get('indel_only_overlap_pairs_post', 0))
                + float(d.get('snv_indel_overlap_pairs_post', 0))
            )
            if pre_total > 0:
                indel_overall = 100.0 * post_total / pre_total
        if indel_overall is not None and indel_overall >= 5:
            times['wasp2rust_indel'] = d['total_s']

    # WASP1
    wasp1_file = first_existing(
        FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp1.json',
        REPO_ROOT / 'benchmarking/star_wasp_comparison/results/wasp1_2025-12-06_20-28-15/benchmark_results.json',
    )
    if wasp1_file is not None and wasp1_file.exists():
        d = json.loads(wasp1_file.read_text())
        times['wasp1'] = d['total_s']

    if len(times) >= 3:
        pipeline_order = ['star_wasp', 'wasp2rust_snv', 'wasp2rust_indel', 'wasp2python', 'wasp1']
        # Use 3-line Rust labels to reduce horizontal crowding in narrow panels
        pipeline_labels = [
            'STAR+\nWASP',
            'WASP2-\nRust\n(SNV)',
            'WASP2-\nRust\n(+INDEL)',
            'WASP2-\nPython',
            'WASP1',
        ]
        colors = ['#56B4E9', C['rust_snv'], C['singlecell'], '#3776AB', C['wasp1']]

        # Filter to only pipelines we have data for
        valid_pipelines = [p for p in pipeline_order if p in times]
        valid_labels = [pipeline_labels[pipeline_order.index(p)] for p in valid_pipelines]
        valid_colors = [colors[pipeline_order.index(p)] for p in valid_pipelines]

        x = np.arange(len(valid_pipelines))
        width = 0.6

        # Total time in minutes
        times_min = np.array([times[p] / 60 for p in valid_pipelines])

        bars = ax.bar(x, times_min, width, color=valid_colors, edgecolor='none')

        # Add time labels on bars
        for i, (bar, time_min) in enumerate(zip(bars, times_min)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   f'{time_min:.1f}m',
                   ha='center', fontsize=4.5, color=C['text'], fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels(valid_labels, fontsize=4.5)
        ax.set_ylabel('Time (min)', fontsize=5)
        ax.set_ylim(bottom=0, top=max(times_min) * 1.15)
        ax.set_title('RNA-seq (HG00731)', fontsize=6, fontweight='bold')

    ax.text(-0.15, 1.05, 'd', transform=ax.transAxes, fontsize=8, fontweight='bold')


def panel_e_rnaseq_filtering(ax):
    """Panel E: RNA-seq retention of variant-overlapping read pairs.

    Definition matches Panel C:
    - "Original" = baseline-filtered variant-overlapping pairs entering WASP remap+filter.
    - "Retained" = those pairs kept after WASP remap filtering.
    Units: read pairs (fragments), stacked by overlap type of the ORIGINAL fragment.
    """
    import json
    from matplotlib.patches import Patch

    type_colors = {
        'snv_only': '#56B4E9',
        'indel_only': '#CC79A7',
        'both': '#009E73',
    }

    def load_split_retention(run_dir, json_path):
        """Return ({pre_by_type}, {post_by_type}) in PAIRS using old-style definitions, or None."""
        if json_path is None or not json_path.exists():
            return None
        d = json.loads(json_path.read_text())
        if 'snv_only_overlap_pairs_pre' not in d:
            return None

        pre = {
            'snv_only': float(d.get('snv_only_overlap_pairs_pre', 0)),
            'indel_only': float(d.get('indel_only_overlap_pairs_pre', 0)),
            'both': float(d.get('snv_indel_overlap_pairs_pre', 0)),
        }
        post = {
            'snv_only': float(d.get('snv_only_overlap_pairs_post', 0)),
            'indel_only': float(d.get('indel_only_overlap_pairs_post', 0)),
            'both': float(d.get('snv_indel_overlap_pairs_post', 0)),
        }

        # Optional safety rescale: ensure pre totals match old-style (original - keep) if present
        if 'original_reads' in d and 'keep_reads' in d:
            old_pre_pairs = (float(d['original_reads']) - float(d['keep_reads'])) / 2.0
            pre_total = sum(pre.values())
            if old_pre_pairs > 0 and pre_total > 0:
                frac_diff = abs(pre_total - old_pre_pairs) / old_pre_pairs
                if frac_diff > 0.02:
                    scale = old_pre_pairs / pre_total
                    pre = {k: v * scale for k, v in pre.items()}

        return pre, post

    results = {}

    # STAR+WASP (old-style split keys backfilled into its JSON)
    star_wasp_file = first_existing(
        FIG1_BENCH_DIR / 'panel_de_rnaseq_star_wasp.json',
        REPO_ROOT / 'benchmarking/star_wasp_comparison/results/star_wasp_2025-12-04_00-58-54/benchmark_results.json',
    )
    if star_wasp_file is not None and star_wasp_file.exists():
        run_dir = star_wasp_file.parent
        r = load_split_retention(run_dir, star_wasp_file)
        if r:
            results['star_wasp'] = {'pre': r[0], 'post': r[1]}

    # WASP2-Python FIXED (latest)
    py_file = first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2python.json')
    if py_file is None:
        py_fixed_files = sorted(
            (REPO_ROOT / 'benchmarking/star_wasp_comparison/results').glob('wasp2python_rnaseq_fixed_*/benchmark_results.json')
        )
        py_file = py_fixed_files[-1] if py_fixed_files else None
    if py_file is not None and py_file.exists():
        run_dir = py_file.parent
        r = load_split_retention(run_dir, py_file)
        if r:
            results['wasp2python'] = {'pre': r[0], 'post': r[1]}

    # WASP2-Rust SNV FAIR (specific SNP-only run)
    rust_snp_fair_file = first_existing(
        FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_snv.json',
        REPO_ROOT / 'benchmarking/star_wasp_comparison/results/wasp2rust_fair_2025-12-10_19-23-02/benchmark_results.json',
    )
    if rust_snp_fair_file is not None and rust_snp_fair_file.exists():
        run_dir = rust_snp_fair_file.parent
        r = load_split_retention(run_dir, rust_snp_fair_file)
        if r:
            results['wasp2rust_snv'] = {'pre': r[0], 'post': r[1]}

    # WASP2-Rust INDEL RNA-seq v2 (latest)
    indel_file = first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_indel.json')
    if indel_file is None:
        rust_indel_files = sorted(
            (REPO_ROOT / 'benchmarking/star_wasp_comparison/results').glob('wasp2rust_indel_rnaseq_v2_*/benchmark_results.json')
        )
        indel_file = rust_indel_files[-1] if rust_indel_files else None
    if indel_file is not None and indel_file.exists():
        run_dir = indel_file.parent
        r = load_split_retention(run_dir, indel_file)
        if r:
            results['wasp2rust_indel'] = {'pre': r[0], 'post': r[1]}

    if len(results) >= 2:
        pipeline_order = []
        pipeline_labels = []
        if 'star_wasp' in results:
            pipeline_order.append('star_wasp')
            pipeline_labels.append('STAR+\nWASP')
        if 'wasp2rust_snv' in results:
            pipeline_order.append('wasp2rust_snv')
            pipeline_labels.append('WASP2-\nRust\n(SNV)')
        if 'wasp2rust_indel' in results:
            pipeline_order.append('wasp2rust_indel')
            pipeline_labels.append('WASP2-\nRust\n(+INDEL)')
        if 'wasp2python' in results:
            pipeline_order.append('wasp2python')
            pipeline_labels.append('WASP2-\nPython')

        x = np.arange(len(pipeline_order))
        width = 0.35

        pre_totals = np.array([sum(results[p]['pre'].values()) for p in pipeline_order], dtype=float)
        post_totals = np.array([sum(results[p]['post'].values()) for p in pipeline_order], dtype=float)
        pass_rates = np.where(pre_totals > 0, 100.0 * post_totals / pre_totals, 0.0)

        bottom_pre = np.zeros(len(pipeline_order))
        bottom_post = np.zeros(len(pipeline_order))
        for key in ['snv_only', 'indel_only', 'both']:
            vals_pre = np.array([results[p]['pre'][key] / 1e6 for p in pipeline_order])
            vals_post = np.array([results[p]['post'][key] / 1e6 for p in pipeline_order])
            ax.bar(
                x - width / 2,
                vals_pre,
                width,
                bottom=bottom_pre,
                color=type_colors[key],
                alpha=0.35,
                edgecolor='none',
            )
            ax.bar(
                x + width / 2,
                vals_post,
                width,
                bottom=bottom_post,
                color=type_colors[key],
                alpha=0.95,
                edgecolor='none',
            )
            bottom_pre += vals_pre
            bottom_post += vals_post

        ax.bar(
            x - width / 2,
            pre_totals / 1e6,
            width,
            bottom=0,
            facecolor='none',
            edgecolor=C['text_light'],
            hatch='///',
            linewidth=0.0,
            zorder=3,
        )

        for i, (post_val, rate) in enumerate(zip(post_totals / 1e6, pass_rates)):
            ax.text(i + width/2, post_val + 0.10, f"{rate:.1f}%",
                    ha='center', fontsize=4.5, color=C['text'], fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels(pipeline_labels, fontsize=4.5)
        ax.set_ylabel('Variant-overlapping\nread pairs (M)', fontsize=5)
        type_handles = [
            Patch(facecolor=type_colors['snv_only'], label='SNV-only'),
            Patch(facecolor=type_colors['indel_only'], label='INDEL-only'),
            Patch(facecolor=type_colors['both'], label='SNV+INDEL'),
        ]
        ax.legend(handles=type_handles, fontsize=3.8, loc='upper right', framealpha=0.95)
        ax.text(
            0.02,
            0.98,
            'Hatched = Original\nSolid = Retained',
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=3.6,
            color=C['text_light'],
        )
        ax.set_ylim(bottom=0, top=max(pre_totals / 1e6) * 1.25)
        ax.set_title('RNA-seq (HG00731)', fontsize=6, fontweight='bold')

    ax.text(-0.15, 1.05, 'e', transform=ax.transAxes, fontsize=8, fontweight='bold')


def generate_figure1():
    """Generate Figure 1 v13 with 5-panel layout (A-E)."""
    setup_style()

    # Figure: 180mm x ~180mm (compact layout with less whitespace)
    fig = plt.figure(figsize=(7.09, 6.5))

    # GridSpec - 2 rows: Panel A (full width), then 4 panels below
    # Reduced height_ratio for Panel A to eliminate whitespace below diagram
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1.2, 1],
                  hspace=0.25, wspace=0.35,
                  left=0.06, right=0.98, top=0.97, bottom=0.06)

    # Panel A: Architecture (top, full width)
    ax_a = fig.add_subplot(gs[0, :])
    panel_a_architecture(ax_a)

    # Panel B: ATAC-seq Benchmark (bottom-left)
    ax_b = fig.add_subplot(gs[1, 0])
    panel_b_benchmark(ax_b)

    # Panel C: ATAC-seq Filtering (bottom-center-left)
    ax_c = fig.add_subplot(gs[1, 1])
    panel_c_filtering(ax_c)

    # Panel D: RNA-seq Benchmark (bottom-center-right)
    ax_d = fig.add_subplot(gs[1, 2])
    panel_d_rnaseq_benchmark(ax_d)

    # Panel E: RNA-seq Filtering (bottom-right)
    ax_e = fig.add_subplot(gs[1, 3])
    panel_e_rnaseq_filtering(ax_e)

    # Save
    out_dir = get_plot_path(1, 'figure1').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    for ext in ['png', 'pdf', 'tiff']:
        plt.savefig(get_plot_path(1, 'figure1', ext),
                   bbox_inches='tight', facecolor='white', dpi=300)
        print(f"Saved: {get_plot_path(1, 'figure1', ext)}")

    return fig


if __name__ == '__main__':
    generate_figure1()
