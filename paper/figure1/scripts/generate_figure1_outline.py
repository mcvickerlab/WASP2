#!/usr/bin/env python3
"""
Figure 1: WASP2 Read Mapping (Paper Outline Version)

Panel A: Schematic of WASP2 pipeline
Panel B: Speed comparison of WASP2 vs. WASP1 remapping for ATAC-seq and RNA-seq (+STAR)
Panel C: Number of reads mapped pre- and post-remapping that overlap (1) SNVs, (2) Indels

Nature Methods compliant: 180mm width, 300dpi, 5-7pt fonts.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
import pandas as pd
import json
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


# Colorblind-safe palette (Paul Tol bright)
C = {
    'comp1': '#4477AA',         # Blue - Read Filtering
    'comp1_bg': '#DCE9F4',
    'comp2': '#228833',         # Green - Counting
    'comp2_bg': '#D9EDDB',
    'comp3': '#CCBB44',         # Yellow/Gold - AI Detection
    'comp3_bg': '#F5F2DC',
    'container_bg': '#F5F8FC',
    'container_border': '#B8C9DC',
    'input': '#555555',
    'input_bg': '#E8E8E8',
    'text': '#222222',
    'text_light': '#555555',
    'arrow': '#666666',
    'white': '#FFFFFF',
    'wasp1': '#888888',
    'wasp2_python': '#3776AB',  # Python blue
    'wasp2_rust': '#E69F00',    # Orange
    'star_wasp': '#56B4E9',     # Sky blue
    'vermillion': '#D55E00',    # Okabe-Ito vermillion (colorblind-safe)
    'snv_only': '#56B4E9',      # Okabe-Ito sky blue
    'indel_only': '#CC79A7',    # Okabe-Ito reddish purple
    'both': '#009E73',          # Okabe-Ito bluish green
}


def setup_style():
    """Nature Methods style."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'axes.titlesize': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.5,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,
        'figure.dpi': 150,
        'savefig.dpi': 300,
    })


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
    box.set_antialiased(True)
    ax.add_patch(box)
    return box


def draw_number_badge(ax, x, y, num, color, radius=0.025):
    """Draw a numbered circle badge."""
    circle = Circle((x, y), radius, facecolor=color, edgecolor='white', lw=1.5, zorder=10)
    ax.add_patch(circle)
    ax.text(x, y, str(num), ha='center', va='center',
            fontsize=7, fontweight='bold', color='white', zorder=11)


def draw_arrow(ax, start, end, color='#666666', lw=1.5):
    """Draw an arrow."""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                               mutation_scale=12, antialiased=True))


def panel_a_schematic(ax):
    """Panel A: WASP2 pipeline schematic.

    Clean 3-component horizontal layout showing the WASP2 workflow:
    Input (BAM/VCF) -> Read Filtering -> Variant Counting -> AI Detection -> Output
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Container
    draw_rounded_box(ax, 0.02, 0.05, 0.96, 0.90, C['container_bg'], C['container_border'], 1.5, alpha=0.5)

    # Title with logo
    if LOGO_PATH.exists():
        logo_img = mpimg.imread(str(LOGO_PATH))
        imagebox = OffsetImage(logo_img, zoom=0.06)
        ab = AnnotationBbox(imagebox, (0.08, 0.88), frameon=False, zorder=5)
        ax.add_artist(ab)
        ax.text(0.14, 0.88, 'WASP2', ha='left', va='center', fontsize=10, fontweight='bold', color=C['comp1'])
    else:
        ax.text(0.08, 0.88, 'WASP2', ha='left', va='center', fontsize=10, fontweight='bold', color=C['comp1'])

    # Input boxes (left side)
    input_x = 0.05
    input_w = 0.10
    input_h = 0.12

    draw_rounded_box(ax, input_x, 0.55, input_w, input_h, C['input_bg'], C['input'], 1.0)
    ax.text(input_x + input_w/2, 0.55 + input_h/2, 'BAM', ha='center', va='center', fontsize=7, fontweight='bold')

    draw_rounded_box(ax, input_x, 0.35, input_w, input_h, C['input_bg'], C['input'], 1.0)
    ax.text(input_x + input_w/2, 0.35 + input_h/2, 'VCF', ha='center', va='center', fontsize=7, fontweight='bold')

    # Arrow from inputs
    draw_arrow(ax, (input_x + input_w + 0.02, 0.50), (0.22, 0.50), C['arrow'], 2.0)

    # Three main components
    comp_y = 0.20
    comp_h = 0.58
    comp_w = 0.20
    gap = 0.06

    comp1_x = 0.24
    comp2_x = comp1_x + comp_w + gap
    comp3_x = comp2_x + comp_w + gap

    # Component 1: Read Filtering
    draw_rounded_box(ax, comp1_x, comp_y, comp_w, comp_h, C['comp1_bg'], C['comp1'], 1.5)
    draw_number_badge(ax, comp1_x + 0.03, comp_y + comp_h - 0.04, 1, C['comp1'])
    ax.text(comp1_x + 0.065, comp_y + comp_h - 0.04, 'Read Filtering', ha='left', va='center',
            fontsize=7, fontweight='bold', color=C['comp1'])

    # Substeps
    substeps1 = ['Identify Variants', 'Swap Alleles', 'Remap Reads', 'Filter Discordant']
    for i, step in enumerate(substeps1):
        sy = comp_y + comp_h - 0.16 - i * 0.11
        draw_rounded_box(ax, comp1_x + 0.015, sy, comp_w - 0.03, 0.08, C['comp1'], None, 0, 0.8)
        ax.text(comp1_x + comp_w/2, sy + 0.04, step, ha='center', va='center', fontsize=5.5, color=C['white'], fontweight='bold')

    ax.text(comp1_x + comp_w/2, comp_y + 0.03, 'SNVs + INDELs', ha='center', va='center', fontsize=5, style='italic', color=C['text_light'])

    # Arrow 1 -> 2
    draw_arrow(ax, (comp1_x + comp_w + 0.015, comp_y + comp_h/2), (comp2_x - 0.015, comp_y + comp_h/2), C['arrow'], 2.0)

    # Component 2: Variant Counting
    draw_rounded_box(ax, comp2_x, comp_y, comp_w, comp_h, C['comp2_bg'], C['comp2'], 1.5)
    draw_number_badge(ax, comp2_x + 0.03, comp_y + comp_h - 0.04, 2, C['comp2'])
    ax.text(comp2_x + 0.065, comp_y + comp_h - 0.04, 'Counting', ha='left', va='center',
            fontsize=7, fontweight='bold', color=C['comp2'])

    substeps2 = ['Het Variants', 'Feature Overlap', 'Allele Counts']
    for i, step in enumerate(substeps2):
        sy = comp_y + comp_h - 0.16 - i * 0.11
        draw_rounded_box(ax, comp2_x + 0.015, sy, comp_w - 0.03, 0.08, C['comp2'], None, 0, 0.8)
        ax.text(comp2_x + comp_w/2, sy + 0.04, step, ha='center', va='center', fontsize=5.5, color=C['white'], fontweight='bold')

    ax.text(comp2_x + comp_w/2, comp_y + 0.03, 'Per-SNV & feature', ha='center', va='center', fontsize=5, style='italic', color=C['text_light'])

    # Arrow 2 -> 3
    draw_arrow(ax, (comp2_x + comp_w + 0.015, comp_y + comp_h/2), (comp3_x - 0.015, comp_y + comp_h/2), C['arrow'], 2.0)

    # Component 3: AI Detection
    draw_rounded_box(ax, comp3_x, comp_y, comp_w, comp_h, C['comp3_bg'], C['comp3'], 1.5)
    draw_number_badge(ax, comp3_x + 0.03, comp_y + comp_h - 0.04, 3, C['comp3'])
    ax.text(comp3_x + 0.065, comp_y + comp_h - 0.04, 'AI Detection', ha='left', va='center',
            fontsize=7, fontweight='bold', color=C['comp3'])

    substeps3 = ['Aggregate', 'Beta-Binomial', 'FDR Correction']
    for i, step in enumerate(substeps3):
        sy = comp_y + comp_h - 0.16 - i * 0.11
        draw_rounded_box(ax, comp3_x + 0.015, sy, comp_w - 0.03, 0.08, C['comp3'], None, 0, 0.8)
        ax.text(comp3_x + comp_w/2, sy + 0.04, step, ha='center', va='center', fontsize=5.5, color=C['white'], fontweight='bold')

    ax.text(comp3_x + comp_w/2, comp_y + 0.03, 'Multi-SNV agg.', ha='center', va='center', fontsize=5, style='italic', color=C['text_light'])

    # Output (right side)
    output_x = 0.90
    draw_arrow(ax, (comp3_x + comp_w + 0.015, comp_y + comp_h/2), (output_x - 0.02, comp_y + comp_h/2), C['arrow'], 2.0)

    draw_rounded_box(ax, output_x - 0.04, comp_y + comp_h/2 - 0.08, 0.10, 0.16, C['white'], C['input'], 1.2)
    ax.text(output_x + 0.01, comp_y + comp_h/2, 'Results', ha='center', va='center', fontsize=7, fontweight='bold')

    ax.text(-0.02, 1.02, 'A', transform=ax.transAxes, fontsize=8, fontweight='bold')


def panel_b_speed_comparison(ax):
    """Panel B: Speed comparison - WASP2 vs WASP1 for ATAC-seq and RNA-seq.

    Grouped bar chart comparing:
    - WASP1
    - WASP2-Python
    - WASP2-Rust (SNV)
    - WASP2-Rust (+INDEL)
    - STAR+WASP (RNA-seq only)

    For both ATAC-seq and RNA-seq datasets.
    """
    # Load benchmark data
    data = {'atac': {}, 'rnaseq': {}}

    # ATAC-seq timing data (from scaling files - use largest size point)
    rust_snp_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp2rust_snv_scaling.tsv',
    )
    rust_indel_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp2rust_indel_scaling.tsv',
    )
    wasp1_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp1_perf_log.tsv',
    )
    wasp2_py_file = first_existing(
        FIG1_BENCH_DIR / 'panel_b_wasp2python_perf_log.tsv',
    )

    # ATAC-seq data
    if rust_snp_file and rust_snp_file.exists():
        df = pd.read_csv(rust_snp_file, sep='\t')
        max_reads = df['n_reads'].max()
        data['atac']['wasp2_rust_snv'] = df[df['n_reads'] == max_reads]['total_s'].mean()

    if rust_indel_file and rust_indel_file.exists():
        df = pd.read_csv(rust_indel_file, sep='\t')
        max_reads = df['n_reads'].max()
        data['atac']['wasp2_rust_indel'] = df[df['n_reads'] == max_reads]['total_s'].mean()

    if wasp1_file and wasp1_file.exists():
        df = pd.read_csv(wasp1_file, sep='\t',
                        names=['n_reads', 'seed', 'total', 'snp2h5', 'intersect', 'remap', 'filter', 'extra'])
        max_reads = df['n_reads'].max()
        data['atac']['wasp1'] = df[df['n_reads'] == max_reads]['total'].mean()

    if wasp2_py_file and wasp2_py_file.exists():
        df = pd.read_csv(wasp2_py_file, sep='\t',
                        names=['n_reads', 'seed', 'total', 'intersect', 'remap', 'filter'])
        max_reads = df['n_reads'].max()
        data['atac']['wasp2_python'] = df[df['n_reads'] == max_reads]['total'].mean()

    # RNA-seq timing data
    rnaseq_files = {
        'star_wasp': first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_star_wasp.json'),
        'wasp2_rust_snv': first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_snv.json'),
        'wasp2_rust_indel': first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_indel.json'),
        'wasp2_python': first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2python.json'),
        'wasp1': first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp1.json'),
    }

    for key, fpath in rnaseq_files.items():
        if fpath and fpath.exists():
            d = json.loads(fpath.read_text())
            time_s = d.get('total_s') or d.get('wall_clock_s', 0)
            if time_s > 0:
                data['rnaseq'][key] = time_s

    # Create grouped bar chart
    methods = ['WASP1', 'WASP2-\nPython', 'WASP2-Rust\n(SNV)', 'WASP2-Rust\n(+INDEL)', 'STAR+\nWASP']
    method_keys = ['wasp1', 'wasp2_python', 'wasp2_rust_snv', 'wasp2_rust_indel', 'star_wasp']
    colors = [C['wasp1'], C['wasp2_python'], C['wasp2_rust'], C['vermillion'], C['star_wasp']]

    x = np.arange(len(methods))
    width = 0.35

    # ATAC-seq bars (left)
    atac_times = []
    for key in method_keys:
        if key == 'star_wasp':
            atac_times.append(0)  # STAR+WASP is RNA-seq only
        else:
            atac_times.append(data['atac'].get(key, 0))

    # RNA-seq bars (right)
    rnaseq_times = []
    for key in method_keys:
        rnaseq_times.append(data['rnaseq'].get(key, 0))

    # Convert to minutes for display
    atac_times_min = np.array(atac_times) / 60
    rnaseq_times_min = np.array(rnaseq_times) / 60

    # Only plot methods with data
    mask_atac = np.array(atac_times) > 0
    mask_rnaseq = np.array(rnaseq_times) > 0

    if any(mask_atac):
        bars1 = ax.bar(x[mask_atac] - width/2, atac_times_min[mask_atac], width,
                       color=[colors[i] for i in range(len(colors)) if mask_atac[i]],
                       alpha=0.7, label='ATAC-seq', edgecolor='black', linewidth=0.5)
        for bar, val in zip(bars1, atac_times_min[mask_atac]):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.1f}',
                       ha='center', va='bottom', fontsize=5, rotation=0)

    if any(mask_rnaseq):
        bars2 = ax.bar(x[mask_rnaseq] + width/2, rnaseq_times_min[mask_rnaseq], width,
                       color=[colors[i] for i in range(len(colors)) if mask_rnaseq[i]],
                       alpha=1.0, label='RNA-seq', edgecolor='black', linewidth=0.5,
                       hatch='///')
        for bar, val in zip(bars2, rnaseq_times_min[mask_rnaseq]):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.1f}',
                       ha='center', va='bottom', fontsize=5, rotation=0)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=6)
    ax.set_ylabel('Time (minutes)', fontsize=7)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=6, loc='upper right',
              facecolor='white', edgecolor='#666666', framealpha=0.9,
              frameon=True, borderpad=0.4)

    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=8, fontweight='bold')
    ax.set_title('Speed Comparison (8 threads)', fontsize=8, fontweight='bold')


def panel_c_read_counts(ax):
    """Panel C: Read counts pre/post remapping for SNVs and INDELs.

    Stacked bar chart showing:
    - Original vs Retained reads
    - Stratified by overlap type (SNV-only, INDEL-only, SNV+INDEL)
    - For both ATAC-seq and RNA-seq
    """
    from matplotlib.patches import Patch

    results = {'atac': {}, 'rnaseq': {}}

    def load_retention(json_path):
        """Load retention data from benchmark JSON."""
        if not json_path or not json_path.exists():
            return None
        d = json.loads(json_path.read_text())

        # New split fields (reads or pairs)
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
            return pre, post
        elif 'snv_only_overlap_pairs_pre' in d:
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
            return pre, post
        return None

    # ATAC-seq retention (from WASP2-Rust INDEL run)
    atac_file = first_existing(FIG1_BENCH_DIR / 'panel_c_atac_wasp2rust_indel.json')
    r = load_retention(atac_file)
    if r:
        results['atac'] = {'pre': r[0], 'post': r[1]}

    # RNA-seq retention (from WASP2-Rust INDEL run)
    rnaseq_file = first_existing(FIG1_BENCH_DIR / 'panel_de_rnaseq_wasp2rust_indel.json')
    r = load_retention(rnaseq_file)
    if r:
        results['rnaseq'] = {'pre': r[0], 'post': r[1]}

    # Plot
    type_colors = {
        'snv_only': C['snv_only'],
        'indel_only': C['indel_only'],
        'both': C['both'],
    }

    datasets = []
    labels = []
    if results.get('atac'):
        datasets.append('atac')
        labels.append('ATAC-seq')
    if results.get('rnaseq'):
        datasets.append('rnaseq')
        labels.append('RNA-seq')

    if not datasets:
        ax.text(0.5, 0.5, 'Read count data needed\n(run benchmarks)',
                transform=ax.transAxes, ha='center', va='center', fontsize=8, color=C['text_light'],
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='#666666',
                          linewidth=0.5, alpha=0.9))
        ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=8, fontweight='bold')
        return

    n_datasets = len(datasets)
    x = np.arange(n_datasets)
    width = 0.35

    # Calculate totals and pass rates
    pre_totals = np.array([sum(results[d]['pre'].values()) for d in datasets], dtype=float)
    post_totals = np.array([sum(results[d]['post'].values()) for d in datasets], dtype=float)
    pass_rates = np.where(pre_totals > 0, 100.0 * post_totals / pre_totals, 0.0)

    # Plot stacked bars
    bottom_pre = np.zeros(n_datasets)
    bottom_post = np.zeros(n_datasets)

    for key in ['snv_only', 'indel_only', 'both']:
        vals_pre = np.array([results[d]['pre'][key] / 1e6 for d in datasets])
        vals_post = np.array([results[d]['post'][key] / 1e6 for d in datasets])

        # Original (hatched)
        ax.bar(x - width/2, vals_pre, width, bottom=bottom_pre,
               color=type_colors[key], alpha=0.35, edgecolor='none')

        # Retained (solid)
        ax.bar(x + width/2, vals_post, width, bottom=bottom_post,
               color=type_colors[key], alpha=0.95, edgecolor='none')

        bottom_pre += vals_pre
        bottom_post += vals_post

    # Hatch overlay for original bars
    ax.bar(x - width/2, pre_totals / 1e6, width, bottom=0,
           facecolor='none', edgecolor=C['text_light'], hatch='///', linewidth=0, zorder=3)

    # Pass rate labels
    for i, (post_val, rate) in enumerate(zip(post_totals / 1e6, pass_rates)):
        ax.text(i + width/2, post_val + 0.1, f'{rate:.1f}%',
                ha='center', fontsize=6, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel('Variant-overlapping\nreads (millions)', fontsize=7)
    ax.set_ylim(bottom=0, top=max(pre_totals / 1e6) * 1.25)

    # Legend
    type_handles = [
        Patch(facecolor=type_colors['snv_only'], label='SNV-only'),
        Patch(facecolor=type_colors['indel_only'], label='INDEL-only'),
        Patch(facecolor=type_colors['both'], label='SNV+INDEL'),
    ]
    ax.legend(handles=type_handles, fontsize=5, loc='upper right',
              facecolor='white', edgecolor='#666666', framealpha=0.9,
              frameon=True, borderpad=0.4)

    ax.text(0.02, 0.98, 'Hatched = Original\nSolid = Retained',
            transform=ax.transAxes, ha='left', va='top', fontsize=5, color=C['text'],
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#666666',
                      linewidth=0.5, alpha=0.9))

    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=8, fontweight='bold')
    ax.set_title('Read Retention', fontsize=8, fontweight='bold')


def generate_figure1():
    """Generate Figure 1 with 3-panel layout matching paper outline."""
    setup_style()

    # Figure: 180mm max width = 7.09 inches
    fig = plt.figure(figsize=(7.0, 5.4))  # 177.8mm - under 180mm limit

    # GridSpec - 2 rows: Panel A (full width), Panels B & C (bottom)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.2, 1],
                  hspace=0.35, wspace=0.35,
                  left=0.06, right=0.98, top=0.95, bottom=0.08)

    # Panel A: Schematic (top, full width)
    ax_a = fig.add_subplot(gs[0, :])
    panel_a_schematic(ax_a)

    # Panel B: Speed comparison (bottom-left)
    ax_b = fig.add_subplot(gs[1, 0])
    panel_b_speed_comparison(ax_b)

    # Panel C: Read counts (bottom-right)
    ax_c = fig.add_subplot(gs[1, 1])
    panel_c_read_counts(ax_c)

    # Save
    out_dir = get_plot_path(1, 'figure1').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    for ext in ['png', 'pdf']:
        outpath = out_dir / f'figure1_outline.{ext}'
        plt.savefig(outpath, bbox_inches='tight', facecolor='white', dpi=300)
        print(f"Saved: {outpath}")

    return fig


if __name__ == '__main__':
    generate_figure1()
