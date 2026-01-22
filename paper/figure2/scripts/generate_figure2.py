#!/usr/bin/env python3
"""
Figure 2: Allele Counting Benchmarks

Panel A: Speed comparison (WASP2-Rust, GATK ASEReadCounter, phASER)
Panel B: Count correlation between tools
Panel C: Bias reduction (original vs remapped BAM)

Updated to support both HG00731 RNA-seq and GM12878 ATAC-seq datasets.
"""
import sys
from pathlib import Path

# Add paper directory to path (works both when run as a script and when imported).
_this_file = globals().get("__file__")
if _this_file:
    _paper_dir = Path(_this_file).resolve().parents[2]  # scripts/ -> figure2/ -> paper/
else:
    # Fallback for environments that import via `python -c` (no __file__).
    # Try common working directories: repo root or within paper/figure2.
    _cwd = Path.cwd().resolve()
    if (_cwd / "paper" / "config.py").exists():
        _paper_dir = _cwd / "paper"
    elif (_cwd / "config.py").exists():
        _paper_dir = _cwd
    else:
        _paper_dir = _cwd
sys.path.insert(0, str(_paper_dir))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
from scipy.stats import pearsonr
from config import COLORS, TOOL_COLORS, PLOT_SETTINGS, get_plot_path, get_data_path
import argparse


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


def panel_a_speed_comparison(ax, dataset='hg00731'):
    """Panel A: Speed comparison across tools."""

    # Load timing data
    timing_file = get_data_path(2, 'timing_results.json')

    if not timing_file.exists():
        # Placeholder
        ax.text(0.5, 0.5, 'Speed Comparison\n(WASP2 vs GATK vs phASER)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.25, f'Data needed:\nRun run_figure2_benchmarks.sh',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')
        ax.set_title('A', fontsize=9, fontweight='bold', loc='left', x=-0.12)
        return

    # Load timing data
    with open(timing_file, 'r') as f:
        timing_data = json.load(f)

    # Get dataset-specific timing
    if dataset == 'hg00731':
        data = timing_data['datasets']['hg00731_rnaseq']['timing']
        dataset_label = 'RNA-seq (HG00731, 56M reads)'
    else:
        data = timing_data['datasets']['gm12878_atacseq']['timing']
        dataset_label = 'ATAC-seq (GM12878, 159M reads)'

    # Extract times
    tools = ['WASP2-Rust', 'GATK', 'phASER']
    times = [
        data['wasp2_rust_s'],
        data['gatk_s'],
        data['phaser_s']
    ]

    # Color scheme
    colors = [TOOL_COLORS['wasp2_rust'], TOOL_COLORS['gatk'], TOOL_COLORS['phaser']]

    x = np.arange(len(tools))
    width = 0.6

    # Plot bars
    bars = ax.bar(x, times, width, color=colors, edgecolor='black', linewidth=0.5)

    # Add time labels on bars
    def format_time(val):
        if val >= 3600:
            return f'{val/3600:.1f}h'
        elif val >= 60:
            return f'{val/60:.1f}m'
        else:
            return f'{val:.0f}s'

    for bar, val in zip(bars, times):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width()/2, val * 1.1, format_time(val),
                    ha='center', va='bottom', fontsize=6)

    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(tools, fontsize=7)
    ax.tick_params(labelsize=6)

    # Add dataset label
    threads = timing_data.get('threads', 8)
    ax.text(0.5, -0.18, f'{dataset_label}, {threads} threads',
            transform=ax.transAxes, ha='center', va='top', fontsize=6, style='italic')

    ax.set_ylabel('Time (seconds)', fontsize=7)
    ax.set_yscale('log')

    # Set y-limits based on data range
    if max(times) > 0:
        ax.set_ylim(min([t for t in times if t > 0]) * 0.5, max(times) * 2)

    # No speedup callouts (keep panel minimal for publication figures).

    ax.set_title('A', fontsize=9, fontweight='bold', loc='left', x=-0.12)


def panel_b_count_comparison(ax, dataset='hg00731'):
    """Panel B: Count correlation between tools (legacy single comparison)."""

    # Load count comparison data
    if dataset == 'hg00731':
        data_file = get_data_path(2, 'hg00731/count_comparison.tsv')
        title_dataset = 'RNA-seq'
    else:
        data_file = get_data_path(2, 'gm12878/count_comparison.tsv')
        title_dataset = 'ATAC-seq'

    if not data_file.exists():
        # Placeholder
        ax.text(0.5, 0.5, 'Count Comparison\n(WASP2 vs GATK)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.25, 'Data needed:\nRun generate_count_comparison.py',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')
        ax.set_title('B', fontsize=9, fontweight='bold', loc='left', x=-0.12)
        return

    df = pd.read_csv(data_file, sep='\t')

    # Filter to variants with counts in both tools
    mask = (df['wasp2_ref'] + df['wasp2_alt'] > 0) & \
           (df['gatk_ref'] + df['gatk_alt'] > 0)

    df_filtered = df[mask]

    if len(df_filtered) == 0:
        ax.text(0.5, 0.5, 'No overlapping counts',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, color='gray')
        ax.set_title('B', fontsize=9, fontweight='bold', loc='left', x=-0.12)
        return

    # Density plot (hexbin) avoids pixelated overplotting in raster exports.
    x_ref = np.log1p(df_filtered["wasp2_ref"])
    y_ref = np.log1p(df_filtered["gatk_ref"])
    x_alt = np.log1p(df_filtered["wasp2_alt"])
    y_alt = np.log1p(df_filtered["gatk_alt"])

    ax.hexbin(x_ref, y_ref, gridsize=60, mincnt=1, bins="log", cmap="Blues", linewidths=0)
    ax.hexbin(x_alt, y_alt, gridsize=60, mincnt=1, bins="log", cmap="OrRd", linewidths=0)

    # Identity line
    max_val = float(
        max(
            x_ref.max(initial=0),
            y_ref.max(initial=0),
            x_alt.max(initial=0),
            y_alt.max(initial=0),
        )
    )
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1.5, zorder=0)

    ax.set_xlabel('WASP2-Rust log(1+count)')
    ax.set_ylabel('GATK log(1+count)')

    # Correlation - combined ref and alt
    all_wasp2 = pd.concat([df_filtered['wasp2_ref'], df_filtered['wasp2_alt']])
    all_gatk = pd.concat([df_filtered['gatk_ref'], df_filtered['gatk_alt']])

    r, p = pearsonr(np.log1p(all_wasp2), np.log1p(all_gatk))

    # Clean panel label - stats in axis labels
    ax.set_title('B', fontsize=9, fontweight='bold', loc='left')
    ax.set_xlabel(f'WASP2 (r={r:.3f}, n={len(df_filtered):,})')


def _scatter_compare_totals(ax, df, x_ref, x_alt, y_ref, y_alt, xlabel, ylabel):
    """Returns (r, n) tuple for use in titles. Returns (None, 0) if no data."""
    x_total = df[x_ref] + df[x_alt]
    y_total = df[y_ref] + df[y_alt]
    mask = (x_total > 0) & (y_total > 0)
    df = df.loc[mask]
    if len(df) == 0:
        ax.text(0.5, 0.5, 'No overlapping counts',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=8, color='gray')
        ax.set_xticks([])
        ax.set_yticks([])
        return None, 0

    x = np.log1p(df[x_ref] + df[x_alt])
    y = np.log1p(df[y_ref] + df[y_alt])
    # Density plot avoids pixelated points and reads better for Nature-style figures.
    ax.hexbin(x, y, gridsize=65, mincnt=1, bins="log", cmap="Blues", linewidths=0)

    max_val = float(max(x.max(), y.max()))
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.4, linewidth=1.5, zorder=0)

    r, _p = pearsonr(x, y)
    # No stats box inside plot - correlation goes in title (set by caller)

    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.tick_params(labelsize=6)
    return r, len(df)


def panel_b_count_comparison_3way(axes, dataset='hg00731'):
    """Panel B: 3-way count comparison between WASP2, GATK, and phASER."""
    if dataset == 'hg00731':
        data_file = get_data_path(2, 'hg00731/count_comparison.tsv')
    else:
        data_file = get_data_path(2, 'gm12878/count_comparison.tsv')

    if not data_file.exists():
        for ax in axes:
            ax.text(0.5, 0.5, 'Data needed:\nRun generate_count_comparison.py',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=7, color='gray')
            ax.set_xticks([])
            ax.set_yticks([])
        axes[0].set_title('B', fontsize=9, fontweight='bold', loc='left', x=-0.12)
        return

    df = pd.read_csv(data_file, sep='\t')

    # 1) WASP2 vs GATK - clean panel label, stats in axis labels
    r1, n1 = _scatter_compare_totals(
        axes[0], df,
        "wasp2_ref", "wasp2_alt",
        "gatk_ref", "gatk_alt",
        "WASP2", "GATK"
    )
    # Clean panel label "B" only - correlation in x-axis label
    axes[0].set_title('B', fontsize=9, fontweight='bold', loc='left', x=-0.12)
    if r1 is not None:
        axes[0].set_xlabel(f'WASP2 (r={r1:.2f})', fontsize=7)

    # 2) WASP2 vs phASER - no panel label, stats in x-axis label
    r2, n2 = _scatter_compare_totals(
        axes[1], df,
        "wasp2_ref", "wasp2_alt",
        "phaser_ref", "phaser_alt",
        "WASP2", "phASER"
    )
    # No title - stats in x-axis label for consistency
    axes[1].set_title('')
    if r2 is not None:
        axes[1].set_xlabel(f'WASP2 (r={r2:.2f})', fontsize=7)

    # 3) GATK vs phASER - no panel label, stats in x-axis label
    r3, n3 = _scatter_compare_totals(
        axes[2], df,
        "gatk_ref", "gatk_alt",
        "phaser_ref", "phaser_alt",
        "GATK", "phASER"
    )
    # No title - stats in x-axis label for consistency
    axes[2].set_title('')
    if r3 is not None:
        axes[2].set_xlabel(f'GATK (r={r3:.2f})', fontsize=7)


def panel_c_bias_reduction(ax, dataset='hg00731'):
    """Panel C: Before/after comparison across methods (original vs WASP2-processed)."""

    if dataset == 'hg00731':
        data_file = get_data_path(2, 'hg00731/before_after_counts.tsv')
    else:
        data_file = get_data_path(2, 'gm12878/before_after_counts.tsv')

    if not data_file.exists():
        ax.text(0.5, 0.55, 'Before/After Count Shift\n(Original vs WASP2-processed)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, style='italic', color='gray')
        ax.text(0.5, 0.28, 'Data needed:\nRun generate_before_after_counts.py',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='gray')
        ax.set_title('C', fontsize=9, fontweight='bold', loc='left', x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]
    df = df[df["method"].isin(methods)]

    # Collect per-method delta distributions
    deltas = [df.loc[df["method"] == m, "delta_ref_ratio"].values for m in methods]
    ns = [len(d) for d in deltas]

    parts = ax.violinplot(
        deltas,
        showmeans=False,
        showmedians=True,
        showextrema=False,
        widths=0.85,
    )
    method_colors = [TOOL_COLORS["wasp2_rust"], TOOL_COLORS["gatk"], TOOL_COLORS["phaser"]]
    for body, c in zip(parts["bodies"], method_colors):
        body.set_facecolor(c)
        body.set_edgecolor("black")
        body.set_alpha(0.75)
        body.set_linewidth(0.5)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.0)

    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)
    ax.set_xticks(np.arange(1, len(methods) + 1))
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_ylabel("Δ ref ratio (processed − original)")
    # Robust y-limits to avoid a few extreme points dominating the scale.
    all_d = np.concatenate([d.astype(float) for d in deltas if len(d)])
    if all_d.size:
        lo, hi = np.quantile(all_d, [0.01, 0.99])
        pad = max(0.01, float((hi - lo) * 0.10))
        ax.set_ylim(float(lo - pad), float(hi + pad))
    else:
        ax.set_ylim(-0.25, 0.25)

    # Annotation: coverage cutoff + n per method if available
    min_total = None
    stats_path = data_file.parent / "before_after_counts_stats.txt"
    if stats_path.exists():
        for line in stats_path.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                except Exception:
                    min_total = None
                break

    cutoff_note = f"total≥{min_total}" if min_total is not None else "coverage-filtered"

    # Report read-loss/retention and mean Δref/Δalt so this panel is clearly a WASP sanity-check.
    lines = []
    for m, n in zip(methods, ns):
        sub = df[df["method"] == m]
        if len(sub) == 0:
            lines.append(f"{m}: n=0")
            continue
        delta_ref_mean = float((sub["filt_ref"] - sub["orig_ref"]).mean())
        delta_alt_mean = float((sub["filt_alt"] - sub["orig_alt"]).mean())
        retention = (sub["filt_total"] / sub["orig_total"]).to_numpy(dtype=float)
        ret_med = float(np.median(retention)) if retention.size else float("nan")
        abs_before = np.abs(sub["orig_ref_ratio"].to_numpy(dtype=float) - 0.5).mean()
        abs_after = np.abs(sub["filt_ref_ratio"].to_numpy(dtype=float) - 0.5).mean()
        red = ((abs_before - abs_after) / abs_before * 100.0) if abs_before > 0 else 0.0
        lines.append(
            f"{m}: n={n:,}  Δref={delta_ref_mean:.2f}  Δalt={delta_alt_mean:.2f}  "
            f"ret(med)={ret_med:.2f}  |bias|↓{red:.1f}%"
        )
    # Clean panel label - cutoff info in axis label
    ax.set_title('C', fontsize=9, fontweight='bold', loc='left')
    ax.set_xlabel(f'Bias reduction ({cutoff_note})')


def _panel_c_retention(ax, df: pd.DataFrame, methods: list[str]) -> None:
    """Per-site retention ratio distribution: log2(filt_total / orig_total) (per tool)."""
    df = df[df["method"].isin(methods)].copy()
    df["retention"] = (df["filt_total"] / df["orig_total"]).astype(float)
    # Log2 transform makes "more removed than added" visually obvious around 0.
    df = df[df["retention"] > 0].copy()
    df["log2_retention"] = np.log2(df["retention"].astype(float))

    vals = [df.loc[df["method"] == m, "log2_retention"].to_numpy(dtype=float) for m in methods]
    ns = [len(v) for v in vals]

    parts = ax.violinplot(
        vals,
        showmeans=False,
        showmedians=True,
        showextrema=False,
        widths=0.85,
    )
    method_colors = [TOOL_COLORS["wasp2_rust"], TOOL_COLORS["gatk"], TOOL_COLORS["phaser"]]
    for body, c in zip(parts["bodies"], method_colors):
        body.set_facecolor(c)
        body.set_edgecolor("black")
        body.set_alpha(0.75)
        body.set_linewidth(0.5)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.0)

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)
    ax.set_xticks(np.arange(1, len(methods) + 1))
    ax.set_xticklabels([])  # bottom panel owns the x labels
    ax.set_ylabel("log2 retention\n(filt/orig)")

    # Robust y-limits
    all_r = np.concatenate([v for v in vals if len(v)])
    if all_r.size:
        lo, hi = np.quantile(all_r, [0.01, 0.99])
        pad = max(0.02, float((hi - lo) * 0.10))
        ax.set_ylim(float(lo - pad), float(hi + pad))
    else:
        ax.set_ylim(-1.5, 1.5)

    # Stats removed from inside plot - key info in y-axis label

    # Inset: fraction of sites decreased/unchanged/increased (based on totals).
    inset = ax.inset_axes([0.08, 0.06, 0.42, 0.32])
    inset.set_title("Sites (Δ total)", fontsize=6.5, pad=2)
    cats = []
    for m in methods:
        sub = df[df["method"] == m]
        if len(sub) == 0:
            cats.append((0.0, 0.0, 0.0))
            continue
        orig = sub["orig_total"].to_numpy()
        filt = sub["filt_total"].to_numpy()
        dec = float((filt < orig).mean() * 100.0)
        same = float((filt == orig).mean() * 100.0)
        inc = float((filt > orig).mean() * 100.0)
        cats.append((dec, same, inc))

    y = np.arange(len(methods))
    dec = np.array([c[0] for c in cats])
    same = np.array([c[1] for c in cats])
    inc = np.array([c[2] for c in cats])
    inset.barh(y, dec, color=COLORS["blue"], alpha=0.55, edgecolor="black", linewidth=0.3, label="decrease")
    inset.barh(y, same, left=dec, color=COLORS["black"], alpha=0.18, edgecolor="black", linewidth=0.3, label="same")
    inset.barh(y, inc, left=dec + same, color=COLORS["vermillion"], alpha=0.55, edgecolor="black", linewidth=0.3, label="increase")
    inset.set_yticks(y)
    inset.set_yticklabels(methods, fontsize=6)
    inset.set_xlim(0, 100)
    inset.set_xticks([0, 50, 100])
    inset.set_xticklabels(["0", "50", "100"], fontsize=6)
    inset.set_xlabel("% sites", fontsize=6, labelpad=1)
    inset.spines["top"].set_visible(False)
    inset.spines["right"].set_visible(False)
    inset.tick_params(axis="both", length=2, pad=1)


def _panel_c_bias_shift(ax, df: pd.DataFrame, methods: list[str], min_total_note: str | None) -> None:
    """Bias shift distribution: Δ ref ratio (processed − original), per tool."""
    df = df[df["method"].isin(methods)].copy()
    deltas = [df.loc[df["method"] == m, "delta_ref_ratio"].to_numpy(dtype=float) for m in methods]
    ns = [len(d) for d in deltas]

    parts = ax.violinplot(
        deltas,
        showmeans=False,
        showmedians=True,
        showextrema=False,
        widths=0.85,
    )
    method_colors = [TOOL_COLORS["wasp2_rust"], TOOL_COLORS["gatk"], TOOL_COLORS["phaser"]]
    for body, c in zip(parts["bodies"], method_colors):
        body.set_facecolor(c)
        body.set_edgecolor("black")
        body.set_alpha(0.75)
        body.set_linewidth(0.5)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.0)

    ax.axhline(0.0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)
    ax.set_xticks(np.arange(1, len(methods) + 1))
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_ylabel("Δ ref ratio\n(filt − orig)")

    all_d = np.concatenate([d for d in deltas if len(d)])
    if all_d.size:
        lo, hi = np.quantile(all_d, [0.01, 0.99])
        pad = max(0.01, float((hi - lo) * 0.10))
        ax.set_ylim(float(lo - pad), float(hi + pad))
    else:
        ax.set_ylim(-0.25, 0.25)

    # No stats box inside plot - min_total info moved to caller's title if needed


def panel_c_two_plots(ax_top, ax_bottom, dataset: str = "hg00731") -> None:
    """Panel C: two-part sanity check per tool (no forced intersection set)."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/before_after_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/before_after_counts.tsv")

    if not data_file.exists():
        ax_top.text(0.5, 0.5, "Data needed:\nRun generate_before_after_counts.py",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]

    # Coverage cutoff note (matches generation script)
    min_total_note = None
    stats_path = data_file.parent / "before_after_counts_stats.txt"
    if stats_path.exists():
        for line in stats_path.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                    min_total_note = f"per-tool site set\nmin_total≥{min_total}"
                except Exception:
                    min_total_note = None
                break

    _panel_c_retention(ax_top, df, methods)
    _panel_c_bias_shift(ax_bottom, df, methods, min_total_note=min_total_note)
    ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)


def panel_c_summary_bars(ax_top, ax_bottom, dataset: str = "hg00731", n_boot: int = 500, seed: int = 0) -> None:
    """
    Panel C (recommended): two simple summary plots (no distribution violins).

    Top: % sites with decreased/unchanged/increased total counts after WASP processing.
    Bottom: mean Δ ref ratio (processed − original) with bootstrap CI.
    """
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/before_after_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/before_after_counts.tsv")

    if not data_file.exists():
        ax_top.text(0.5, 0.5, "Data needed:\nRun generate_before_after_counts.py",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    # PI intent: Panel C is a WASP filtering sanity check; include filters (WASP2/WASP1),
    # not the counting tools (GATK/phASER) which can confuse interpretation.
    methods = ["WASP2"]
    if "WASP1" in set(df["method"].astype(str)):
        methods.append("WASP1")
    df = df[df["method"].isin(methods)].copy()

    # Coverage cutoff note (matches generation script)
    min_total_note = None
    stats_path = data_file.parent / "before_after_counts_stats.txt"
    if stats_path.exists():
        for line in stats_path.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                    min_total_note = f"per-tool site set\nmin_total≥{min_total}"
                except Exception:
                    min_total_note = None
                break

    # --- Top: site-level total change categories ---
    cats = []
    ns = []
    for m in methods:
        sub = df[df["method"] == m]
        n = int(len(sub))
        ns.append(n)
        if n == 0:
            cats.append((0.0, 0.0, 0.0))
            continue
        orig = sub["orig_total"].to_numpy()
        filt = sub["filt_total"].to_numpy()
        dec = float((filt < orig).mean() * 100.0)
        same = float((filt == orig).mean() * 100.0)
        inc = float((filt > orig).mean() * 100.0)
        cats.append((dec, same, inc))

    y = np.arange(len(methods))
    dec = np.array([c[0] for c in cats])
    same = np.array([c[1] for c in cats])
    inc = np.array([c[2] for c in cats])

    ax_top.barh(y, dec, color=COLORS["blue"], alpha=0.55, edgecolor="black", linewidth=0.4, label="decrease")
    ax_top.barh(y, same, left=dec, color=COLORS["black"], alpha=0.18, edgecolor="black", linewidth=0.4, label="same")
    ax_top.barh(y, inc, left=dec + same, color=COLORS["vermillion"], alpha=0.55, edgecolor="black", linewidth=0.4, label="increase")
    ax_top.set_yticks(y)
    ax_top.set_yticklabels(methods, fontsize=8)
    ax_top.set_xlim(0, 100)
    ax_top.set_xlabel("% SNPs (Δ total counts)", fontsize=8)
    ax_top.legend(fontsize=6.5, framealpha=0.95, loc="lower right")
    ax_top.tick_params(labelsize=7)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    for yi, n in zip(y, ns):
        ax_top.text(0.01, yi, f"n={n:,}", transform=ax_top.get_yaxis_transform(), ha="left", va="center", fontsize=7)

    # No stats box inside plot - min_total info in main figure title if needed

    # --- Bottom: histogram of Δ ref ratio (distribution), per PI feedback ---
    bins = np.linspace(-0.5, 0.5, 61)
    for m in methods:
        sub = df[df["method"] == m]
        x = sub["delta_ref_ratio"].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        color = TOOL_COLORS["wasp2_rust"] if m == "WASP2" else TOOL_COLORS["wasp1"]
        ax_bottom.hist(
            x,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.2,
            color=color,
            label=m,
        )
    ax_bottom.axvline(0.0, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax_bottom.set_xlabel("Δ ref ratio (filt − orig)", fontsize=8)
    ax_bottom.set_ylabel("Density", fontsize=8)
    ax_bottom.tick_params(labelsize=7)
    ax_bottom.legend(fontsize=6.5, framealpha=0.95, loc="upper right")
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)


def panel_c_ref_ratio_before_after(ax_top, ax_bottom, dataset: str = "hg00731") -> None:
    """Alternative Panel C bottom: ref-ratio histograms before vs after filtering (n>=min_total)."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/before_after_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/before_after_counts.tsv")

    if not data_file.exists():
        ax_top.text(0.5, 0.5, "Data needed:\nRun generate_before_after_counts.py",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    # Coverage cutoff note (min_total from generation script)
    min_total_note = None
    stats_path = data_file.parent / "before_after_counts_stats.txt"
    if stats_path.exists():
        for line in stats_path.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                    min_total_note = f"n≥{min_total} pre & post"
                except Exception:
                    min_total_note = None
                break
    methods = ["WASP2"]
    if "WASP1" in set(df["method"].astype(str)):
        methods.append("WASP1")
    df = df[df["method"].isin(methods)].copy()

    # Top: keep the % decrease/same/increase panel, but only for these methods.
    panel_c_summary_bars(ax_top, ax_bottom, dataset=dataset)

    # Overwrite bottom with before/after ref-ratio histograms (instead of Δ histogram).
    ax_bottom.cla()
    bins = np.linspace(0.0, 1.0, 51)
    for m in methods:
        sub = df[df["method"] == m]
        orig = sub["orig_ref_ratio"].to_numpy(dtype=float)
        filt = sub["filt_ref_ratio"].to_numpy(dtype=float)
        orig = orig[np.isfinite(orig)]
        filt = filt[np.isfinite(filt)]
        color = TOOL_COLORS["wasp2_rust"] if m == "WASP2" else TOOL_COLORS["wasp1"]
        ax_bottom.hist(orig, bins=bins, density=True, histtype="step", linewidth=1.2, color=color, alpha=0.5, label=f"{m} orig")
        ax_bottom.hist(filt, bins=bins, density=True, histtype="step", linewidth=1.2, color=color, alpha=1.0, label=f"{m} filt")

    ax_bottom.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax_bottom.set_xlabel("Ref ratio (ref / (ref+alt))", fontsize=8)
    ax_bottom.set_ylabel("Density", fontsize=8)
    ax_bottom.tick_params(labelsize=7)
    ax_bottom.legend(fontsize=6.0, framealpha=0.95, loc="upper right", ncol=2)
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)

    # No stats box inside plot - info in title if needed
    ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)


def panel_c_histograms(ax_top, ax_bottom, dataset: str = "hg00731", min_reads: int = 10) -> None:
    """
    Panel C (per user feedback): clear histograms for WASP2 bias reduction.

    Top: Histogram of ref_ratio BEFORE vs AFTER WASP2 filtering (overlaid).
         Shows the shift toward 0.5 (balanced allelic expression).

    Bottom: Histogram of delta ref_ratio (filtered - original).
            Shows how much each site moved toward/away from balance.

    Only shows WASP2 (the filtering method), not GATK/phASER (counting tools).
    Filters to SNPs with total reads >= min_reads.
    """
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/before_after_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/before_after_counts.tsv")

    if not data_file.exists():
        ax_top.text(0.5, 0.5, "Data needed:\nRun generate_before_after_counts.py",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")

    # Only use WASP2 (the filtering tool) - GATK/phASER are just counting tools
    # and don't perform any filtering, so "before/after" doesn't make sense for them
    df = df[df["method"] == "WASP2"].copy()

    # Apply minimum read filter (user requested n>=10)
    df = df[(df["orig_total"] >= min_reads) & (df["filt_total"] >= min_reads)].copy()

    if len(df) == 0:
        ax_top.text(0.5, 0.5, f"No sites with ≥{min_reads} reads\nbefore AND after filtering",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=8, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    n_sites = len(df)

    # Get ref ratios
    orig_ratio = df["orig_ref_ratio"].to_numpy(dtype=float)
    filt_ratio = df["filt_ref_ratio"].to_numpy(dtype=float)
    delta_ratio = df["delta_ref_ratio"].to_numpy(dtype=float)

    # Remove NaN values
    valid = np.isfinite(orig_ratio) & np.isfinite(filt_ratio) & np.isfinite(delta_ratio)
    orig_ratio = orig_ratio[valid]
    filt_ratio = filt_ratio[valid]
    delta_ratio = delta_ratio[valid]

    # === TOP PANEL: Before vs After ref_ratio histograms ===
    bins_ratio = np.linspace(0.0, 1.0, 41)  # 40 bins from 0 to 1

    # Original (before filtering) - lighter color
    ax_top.hist(orig_ratio, bins=bins_ratio, density=True, histtype="stepfilled",
                alpha=0.4, color=COLORS["black"], edgecolor=COLORS["black"],
                linewidth=1.0, label="Before WASP2")

    # Filtered (after WASP2) - darker color
    ax_top.hist(filt_ratio, bins=bins_ratio, density=True, histtype="stepfilled",
                alpha=0.6, color=TOOL_COLORS["wasp2_rust"], edgecolor=TOOL_COLORS["wasp2_rust"],
                linewidth=1.0, label="After WASP2")

    # Reference line at 0.5 (balanced)
    ax_top.axvline(0.5, color="black", linestyle="--", linewidth=1.2, alpha=0.7)

    ax_top.set_xlabel("Ref ratio (ref / total)", fontsize=7)
    ax_top.set_ylabel("Density", fontsize=7)
    ax_top.tick_params(labelsize=6)
    ax_top.legend(fontsize=6, framealpha=0.95, loc="upper left")
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    # Calculate bias reduction statistics
    bias_before = np.abs(orig_ratio - 0.5).mean()
    bias_after = np.abs(filt_ratio - 0.5).mean()
    bias_reduction_pct = ((bias_before - bias_after) / bias_before * 100) if bias_before > 0 else 0

    # Clean panel label "C" only - stats in annotation box (bottom-right to avoid legend overlap)
    ax_top.set_title('C', fontsize=9, fontweight='bold', loc='left', x=-0.12)
    # Add stats as text annotation in bottom-right (legend is in upper-left)
    ax_top.text(0.98, 0.04, f'n={len(orig_ratio):,}, bias↓{bias_reduction_pct:.0f}%',
                transform=ax_top.transAxes, fontsize=6, ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='none'))

    # === BOTTOM PANEL: Delta ref_ratio histogram ===
    # Symmetric bins around 0
    max_delta = min(0.5, np.quantile(np.abs(delta_ratio), 0.99) * 1.5)
    bins_delta = np.linspace(-max_delta, max_delta, 51)

    ax_bottom.hist(delta_ratio, bins=bins_delta, density=True, histtype="stepfilled",
                   alpha=0.6, color=TOOL_COLORS["wasp2_rust"], edgecolor="black",
                   linewidth=0.5)

    # Reference line at 0 (no change)
    ax_bottom.axvline(0, color="black", linestyle="--", linewidth=1.2, alpha=0.7)

    ax_bottom.set_xlabel(f"Δ ref ratio (n≥{min_reads})", fontsize=7)
    ax_bottom.set_ylabel("Density", fontsize=7)
    ax_bottom.tick_params(labelsize=6)
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)

    # Stats annotation in corner
    mean_delta = np.mean(delta_ratio)
    ax_bottom.text(0.98, 0.98, f'mean Δ={mean_delta:.3f}',
                   transform=ax_bottom.transAxes, fontsize=6, ha='right', va='top',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='none'))

def save_panel_c_raw_counts(dataset: str = "hg00731") -> tuple[Path, Path] | None:
    """Supplementary view: raw ref/alt counts pre vs post (log10 scale), per tool."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/before_after_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/before_after_counts.tsv")

    if not data_file.exists():
        return None

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]
    df = df[df["method"].isin(methods)].copy()

    setup_style()
    # Use constrained_layout for automatic spacing
    fig, axes = plt.subplots(nrows=len(methods), ncols=2, figsize=(6.6, 6.2),
                             sharex=True, sharey=True, layout="constrained")
    if len(methods) == 1:
        axes = np.array([axes])

    bins = np.linspace(0, 3.0, 60)  # log10(1+count) up to ~1000
    for i, m in enumerate(methods):
        sub = df[df["method"] == m]
        for j, allele in enumerate(["ref", "alt"]):
            ax = axes[i, j]
            orig = np.log10(1.0 + sub[f"orig_{allele}"].to_numpy(dtype=float))
            filt = np.log10(1.0 + sub[f"filt_{allele}"].to_numpy(dtype=float))
            ax.hist(orig, bins=bins, density=True, histtype="step", linewidth=1.1, color=COLORS["black"], alpha=0.55, label="original")
            ax.hist(filt, bins=bins, density=True, histtype="step", linewidth=1.1, color=COLORS["blue"], label="WASP2-processed")
            if i == 0:
                ax.set_title(f"{allele.upper()} counts", fontsize=9)
            if j == 0:
                ax.set_ylabel(m, fontsize=9)
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if i == 0 and j == 1:
                ax.legend(fontsize=7, framealpha=0.95)

    for ax in axes[-1, :]:
        ax.set_xlabel("log10(1 + count)", fontsize=8)

    fig.suptitle("Figure 2C (supplement): raw counts", fontsize=10, fontweight="bold")
    # constrained_layout handles spacing automatically

    out_dir = get_plot_path(2, "figure2").parent
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "figure2_c_raw_counts.png"
    pdf_path = out_dir / "figure2_c_raw_counts.pdf"
    fig.savefig(png_path, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png_path, pdf_path


def panel_c_delta_counts(ax, dataset="hg00731"):
    """Panel C: Per-site Δ counts (processed − original) for ref and alt, per method."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/delta_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/delta_counts.tsv")

    if not data_file.exists():
        ax.text(
            0.5,
            0.55,
            "Data needed:\nRun generate_delta_counts.py",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7,
            color="gray",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]
    df = df[df["method"].isin(methods)].copy()
    df["method"] = pd.Categorical(df["method"], categories=methods, ordered=True)

    # Prepare boxplot data
    ref_data = [df.loc[df["method"] == m, "delta_ref"].to_numpy() for m in methods]
    alt_data = [df.loc[df["method"] == m, "delta_alt"].to_numpy() for m in methods]
    ns = [len(x) for x in ref_data]

    positions_ref = np.arange(len(methods)) * 2.0
    positions_alt = positions_ref + 0.7

    def draw_box(data, positions, color, label):
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="black", linewidth=1.0),
            boxprops=dict(edgecolor="black", linewidth=0.5),
            whiskerprops=dict(color="black", linewidth=0.5),
            capprops=dict(color="black", linewidth=0.5),
        )
        for b in bp["boxes"]:
            b.set_facecolor(color)
            b.set_alpha(0.75)
        return bp

    draw_box(ref_data, positions_ref, COLORS["blue"], "Ref")
    draw_box(alt_data, positions_alt, COLORS["vermillion"], "Alt")

    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)
    ax.set_ylabel("Δ allele counts per SNP\n(processed − original)", fontsize=8)
    # Many sites have Δ=0 with a small tail of larger changes; symlog shows both clearly.
    ax.set_yscale("symlog", linthresh=2.0, linscale=1.0)

    # X ticks at method centers
    centers = (positions_ref + positions_alt) / 2.0
    ax.set_xticks(centers)
    ax.set_xticklabels(methods, fontsize=8)

    # Legend and n annotation
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color=COLORS["blue"], lw=6, alpha=0.75, label="Ref"),
            plt.Line2D([0], [0], color=COLORS["vermillion"], lw=6, alpha=0.75, label="Alt"),
        ],
        fontsize=6,
        loc="upper left",
        framealpha=0.95,
    )
    for xc, n in zip(centers, ns):
        ax.text(xc, 0.03, f"n={n:,}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)

    # Use robust y-limits from per-method percentiles if available
    stats_path = data_file.parent / "delta_counts_stats.txt"
    if stats_path.exists():
        lines = stats_path.read_text().splitlines()
        # Use p01/p99 across ref+alt to set an informative plot range while avoiding extreme tails.
        p01_ref = p99_ref = p01_alt = p99_alt = None
        for line in lines:
            if line.startswith("GATK_delta_ref_p01:"):
                p01_ref = float(line.split(":", 1)[1].strip())
            if line.startswith("GATK_delta_ref_p99:"):
                p99_ref = float(line.split(":", 1)[1].strip())
            if line.startswith("GATK_delta_alt_p01:"):
                p01_alt = float(line.split(":", 1)[1].strip())
            if line.startswith("GATK_delta_alt_p99:"):
                p99_alt = float(line.split(":", 1)[1].strip())

        if None not in (p01_ref, p99_ref, p01_alt, p99_alt):
            lo = min(p01_ref, p01_alt)
            hi = max(p99_ref, p99_alt)
            # Expand slightly and snap to a nice round bound.
            lo = float(np.floor(lo * 1.2 / 5.0) * 5.0)
            hi = float(np.ceil(hi * 1.2 / 5.0) * 5.0)
            # Ensure some visible range around zero.
            if hi < 5:
                hi = 5.0
            if lo > -5:
                lo = -5.0
            ax.set_ylim(lo, hi)

    ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)


def panel_c_delta_mean_ci(ax, dataset="hg00731", n_boot: int = 200, seed: int = 0):
    """Panel C style 1: mean Δref/Δalt per SNP with bootstrap CI."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/delta_counts.tsv")
        stats_file = get_data_path(2, "hg00731/delta_counts_stats.txt")
    else:
        data_file = get_data_path(2, "gm12878/delta_counts.tsv")
        stats_file = get_data_path(2, "gm12878/delta_counts_stats.txt")

    if not data_file.exists():
        ax.text(0.5, 0.55, "Data needed:\nRun generate_delta_counts.py",
                transform=ax.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]
    df = df[df["method"].isin(methods)].copy()
    rng = np.random.default_rng(seed)

    def boot_ci(x: np.ndarray) -> tuple[float, float, float]:
        x = x.astype(float)
        n = len(x)
        if n == 0:
            return 0.0, 0.0, 0.0
        means = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, size=n)
            means.append(float(x[idx].mean()))
        means = np.array(means)
        return float(x.mean()), float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))

    rows = []
    for m in methods:
        d = df[df["method"] == m]
        mean_ref, lo_ref, hi_ref = boot_ci(d["delta_ref"].to_numpy())
        mean_alt, lo_alt, hi_alt = boot_ci(d["delta_alt"].to_numpy())
        pct_ref_nz = float((d["delta_ref"] != 0).mean() * 100) if len(d) else 0.0
        pct_alt_nz = float((d["delta_alt"] != 0).mean() * 100) if len(d) else 0.0
        rows.append((m, mean_ref, lo_ref, hi_ref, mean_alt, lo_alt, hi_alt, pct_ref_nz, pct_alt_nz, int(len(d))))

    x = np.arange(len(methods))
    width = 0.28

    mean_ref = np.array([r[1] for r in rows])
    lo_ref = np.array([r[2] for r in rows])
    hi_ref = np.array([r[3] for r in rows])
    mean_alt = np.array([r[4] for r in rows])
    lo_alt = np.array([r[5] for r in rows])
    hi_alt = np.array([r[6] for r in rows])

    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)

    ax.bar(x - width / 2, mean_ref, width, color=COLORS["blue"], edgecolor="black", linewidth=0.4, label="Ref")
    ax.errorbar(x - width / 2, mean_ref, yerr=[mean_ref - lo_ref, hi_ref - mean_ref],
                fmt="none", ecolor="black", elinewidth=0.8, capsize=2, capthick=0.8)

    ax.bar(x + width / 2, mean_alt, width, color=COLORS["vermillion"], edgecolor="black", linewidth=0.4, label="Alt")
    ax.errorbar(x + width / 2, mean_alt, yerr=[mean_alt - lo_alt, hi_alt - mean_alt],
                fmt="none", ecolor="black", elinewidth=0.8, capsize=2, capthick=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_ylabel("Mean Δ counts per SNP\n(processed − original)", fontsize=8)
    ax.legend(fontsize=6, loc="upper left", framealpha=0.95)

    # annotate n and nonzero fractions
    for xi, r in zip(x, rows):
        _m, _mr, _lr, _hr, _ma, _la, _ha, pr, pa, n = r
        ax.text(xi, 0.03, f"n={n:,}\nΔref≠0 {pr:.0f}%\nΔalt≠0 {pa:.0f}%",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5)

    min_total = None
    if stats_file.exists():
        for line in stats_file.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                except Exception:
                    min_total = None
                break
    # Clean panel label only - min_total info in text annotation
    ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
    if min_total is not None:
        ax.text(0.98, 0.98, f'min_total≥{min_total}', transform=ax.transAxes,
                fontsize=6, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))


def panel_c_delta_skew(ax, dataset="hg00731"):
    """Panel C style 2: distribution of Δskew = Δref − Δalt per SNP (count-space bias proxy)."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/delta_counts.tsv")
        stats_file = get_data_path(2, "hg00731/delta_counts_stats.txt")
    else:
        data_file = get_data_path(2, "gm12878/delta_counts.tsv")
        stats_file = get_data_path(2, "gm12878/delta_counts_stats.txt")

    if not data_file.exists():
        ax.text(0.5, 0.55, "Data needed:\nRun generate_delta_counts.py",
                transform=ax.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]
    df = df[df["method"].isin(methods)].copy()
    df["delta_skew"] = df["delta_ref"] - df["delta_alt"]

    data = [df.loc[df["method"] == m, "delta_skew"].to_numpy() for m in methods]
    ns = [len(d) for d in data]

    ax.boxplot(
        data,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.0),
        boxprops=dict(edgecolor="black", linewidth=0.5, facecolor=COLORS["blue"], alpha=0.6),
        whiskerprops=dict(color="black", linewidth=0.5),
        capprops=dict(color="black", linewidth=0.5),
    )

    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)
    ax.set_xticks(np.arange(1, len(methods) + 1))
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_ylabel("Δskew per SNP\n(Δref − Δalt)", fontsize=8)
    ax.set_yscale("symlog", linthresh=2.0, linscale=1.0)

    for i, n in enumerate(ns, start=1):
        ax.text(i, 0.03, f"n={n:,}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)

    min_total = None
    if stats_file.exists():
        for line in stats_file.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                except Exception:
                    min_total = None
                break
    # Clean panel label only - min_total info in text annotation
    ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
    if min_total is not None:
        ax.text(0.98, 0.98, f'min_total≥{min_total}', transform=ax.transAxes,
                fontsize=6, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))


def panel_c_delta_counts_affected_only(ax, dataset="hg00731"):
    """Panel C style 3: same as Δcounts boxplots, but restricted to sites with any change."""
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/delta_counts.tsv")
        stats_file = get_data_path(2, "hg00731/delta_counts_stats.txt")
    else:
        data_file = get_data_path(2, "gm12878/delta_counts.tsv")
        stats_file = get_data_path(2, "gm12878/delta_counts_stats.txt")

    if not data_file.exists():
        ax.text(0.5, 0.55, "Data needed:\nRun generate_delta_counts.py",
                transform=ax.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")
    methods = ["WASP2", "GATK", "phASER"]
    df = df[df["method"].isin(methods)].copy()

    # Identify affected SNPs in the intersection set: any method has any change.
    site_any = (
        df.assign(any_change=(df["delta_ref"] != 0) | (df["delta_alt"] != 0))
        .groupby(["chrom", "pos"], as_index=False)["any_change"]
        .any()
    )
    df = df.merge(site_any, on=["chrom", "pos"], how="inner")
    df = df[df["any_change"]].copy()

    # Reuse the base delta-counts plot, but on this filtered df.
    methods = ["WASP2", "GATK", "phASER"]
    ref_data = [df.loc[df["method"] == m, "delta_ref"].to_numpy() for m in methods]
    alt_data = [df.loc[df["method"] == m, "delta_alt"].to_numpy() for m in methods]
    ns = [len(x) for x in ref_data]

    positions_ref = np.arange(len(methods)) * 2.0
    positions_alt = positions_ref + 0.7

    def draw_box(data, positions, color):
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="black", linewidth=1.0),
            boxprops=dict(edgecolor="black", linewidth=0.5),
            whiskerprops=dict(color="black", linewidth=0.5),
            capprops=dict(color="black", linewidth=0.5),
        )
        for b in bp["boxes"]:
            b.set_facecolor(color)
            b.set_alpha(0.75)

    draw_box(ref_data, positions_ref, COLORS["blue"])
    draw_box(alt_data, positions_alt, COLORS["vermillion"])

    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6, zorder=0)
    ax.set_ylabel("Δ counts per SNP\n(processed − original)\n(affected sites only)", fontsize=8)
    ax.set_yscale("symlog", linthresh=2.0, linscale=1.0)

    centers = (positions_ref + positions_alt) / 2.0
    ax.set_xticks(centers)
    ax.set_xticklabels(methods, fontsize=8)

    ax.legend(
        handles=[
            plt.Line2D([0], [0], color=COLORS["blue"], lw=6, alpha=0.75, label="Ref"),
            plt.Line2D([0], [0], color=COLORS["vermillion"], lw=6, alpha=0.75, label="Alt"),
        ],
        fontsize=6,
        loc="upper left",
        framealpha=0.95,
    )
    for xc, n in zip(centers, ns):
        ax.text(xc, 0.03, f"n={n:,}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)

    min_total = None
    if stats_file.exists():
        for line in stats_file.read_text().splitlines():
            if line.startswith("min_total:"):
                try:
                    min_total = int(line.split(":", 1)[1].strip())
                except Exception:
                    min_total = None
                break
    # Clean panel label only - min_total info in text annotation
    ax.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
    if min_total is not None:
        ax.text(0.98, 0.98, f'min_total≥{min_total}', transform=ax.transAxes,
                fontsize=6, ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))


def panel_c_scatter_stratified(ax_top, ax_bottom, dataset: str = "hg00731", min_reads: int = 10) -> None:
    """
    Panel C: Scatter plot + stratified bias reduction bars.

    Top panel: Hexbin scatter showing ref_ratio before vs after WASP2 filtering.
    Bottom panel: Stratified bias reduction bars showing mean |bias| by initial bias level.

    Parameters:
    - ax_top: matplotlib axis for scatter plot
    - ax_bottom: matplotlib axis for stratified bars
    - dataset: 'hg00731' or 'gm12878'
    - min_reads: minimum total reads required (before AND after filtering)
    """
    if dataset == "hg00731":
        data_file = get_data_path(2, "hg00731/before_after_counts.tsv")
    else:
        data_file = get_data_path(2, "gm12878/before_after_counts.tsv")

    if not data_file.exists():
        ax_top.text(0.5, 0.5, "Data needed:\nRun generate_before_after_counts.py",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=7, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    df = pd.read_csv(data_file, sep="\t")

    # Filter to WASP2 only and sites with sufficient coverage
    df = df[df["method"] == "WASP2"].copy()
    df = df[(df["orig_total"] >= min_reads) & (df["filt_total"] >= min_reads)].copy()

    if len(df) == 0:
        ax_top.text(0.5, 0.5, f"No sites with ≥{min_reads} reads\nbefore AND after filtering",
                    transform=ax_top.transAxes, ha="center", va="center", fontsize=8, color="gray")
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_bottom.set_xticks([])
        ax_bottom.set_yticks([])
        ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
        return

    orig_ratio = df["orig_ref_ratio"].to_numpy(dtype=float)
    filt_ratio = df["filt_ref_ratio"].to_numpy(dtype=float)

    # Remove NaN/inf
    valid = np.isfinite(orig_ratio) & np.isfinite(filt_ratio)
    orig_ratio = orig_ratio[valid]
    filt_ratio = filt_ratio[valid]

    n_sites = len(orig_ratio)

    # === TOP PANEL: Scatter plot (hexbin) ===
    hb = ax_top.hexbin(orig_ratio, filt_ratio, gridsize=50, mincnt=1, bins="log",
                       cmap="Blues", linewidths=0, zorder=1)

    # Identity line (y=x, no change)
    ax_top.plot([0, 1], [0, 1], 'k--', linewidth=1.5, alpha=0.7, label='no change', zorder=2)

    # Reference lines at 0.5 (balanced)
    ax_top.axvline(0.5, color='gray', linestyle=':', linewidth=1.0, alpha=0.5, zorder=0)
    ax_top.axhline(0.5, color='gray', linestyle=':', linewidth=1.0, alpha=0.5, zorder=0)

    ax_top.set_xlabel("Ref ratio before filtering", fontsize=8)
    ax_top.set_ylabel("Ref ratio after filtering", fontsize=8)
    ax_top.set_xlim(0, 1)
    ax_top.set_ylim(0, 1)
    ax_top.tick_params(labelsize=7)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    # Correlation
    r, _p = pearsonr(orig_ratio, filt_ratio)

    # Clean panel label only - stats in text annotation
    ax_top.set_title("C", fontsize=9, fontweight="bold", loc="left", x=-0.12)
    ax_top.text(0.98, 0.02, f'n={n_sites:,}\nr={r:.3f}',
                transform=ax_top.transAxes, fontsize=6, ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='none'))

    # === BOTTOM PANEL: Stratified bias reduction bars ===
    # Calculate absolute bias
    bias_before = np.abs(orig_ratio - 0.5)
    bias_after = np.abs(filt_ratio - 0.5)

    # Define stratification bins based on INITIAL bias
    bins = [0, 0.05, 0.10, 0.20, 0.50]
    bin_labels = ["[0, 0.05)", "[0.05, 0.10)", "[0.10, 0.20)", "[0.20, 0.50]"]

    # Assign each site to a bias category
    bin_indices = np.digitize(bias_before, bins[1:])

    # Calculate statistics for each category
    categories = []
    mean_before = []
    mean_after = []
    err_before = []
    err_after = []
    n_per_cat = []
    pct_reduction = []

    for i, label in enumerate(bin_labels):
        mask = (bin_indices == i)
        n_cat = mask.sum()
        n_per_cat.append(n_cat)

        if n_cat == 0:
            mean_before.append(0)
            mean_after.append(0)
            err_before.append(0)
            err_after.append(0)
            pct_reduction.append(0)
            categories.append(label)
            continue

        before_vals = bias_before[mask]
        after_vals = bias_after[mask]

        # Mean
        m_before = float(before_vals.mean())
        m_after = float(after_vals.mean())
        mean_before.append(m_before)
        mean_after.append(m_after)

        # SEM for error bars
        sem_before = float(before_vals.std() / np.sqrt(n_cat)) if n_cat > 1 else 0
        sem_after = float(after_vals.std() / np.sqrt(n_cat)) if n_cat > 1 else 0
        err_before.append(sem_before)
        err_after.append(sem_after)

        # Percent reduction
        if m_before > 0:
            reduction = ((m_before - m_after) / m_before * 100)
        else:
            reduction = 0
        pct_reduction.append(reduction)

        categories.append(label)

    x = np.arange(len(categories))
    width = 0.35

    # Before bars (lighter)
    bars1 = ax_bottom.bar(x - width/2, mean_before, width,
                          color=COLORS["black"], alpha=0.4,
                          edgecolor="black", linewidth=0.4,
                          label="Before WASP2")
    ax_bottom.errorbar(x - width/2, mean_before, yerr=err_before,
                       fmt="none", ecolor="black", elinewidth=0.6, capsize=2, capthick=0.6)

    # After bars (darker, WASP2 color)
    bars2 = ax_bottom.bar(x + width/2, mean_after, width,
                          color=TOOL_COLORS["wasp2_rust"], alpha=0.7,
                          edgecolor="black", linewidth=0.4,
                          label="After WASP2")
    ax_bottom.errorbar(x + width/2, mean_after, yerr=err_after,
                       fmt="none", ecolor="black", elinewidth=0.6, capsize=2, capthick=0.6)

    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(categories, fontsize=7, rotation=0)
    # Put filter info in x-axis label instead of inside plot
    ax_bottom.set_xlabel(f"Initial |bias| category (sites with >={min_reads} reads)", fontsize=8)
    ax_bottom.set_ylabel("Mean |ref_ratio - 0.5|", fontsize=8)
    ax_bottom.tick_params(labelsize=7)
    ax_bottom.legend(fontsize=6, framealpha=0.95, loc="upper right")
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)


def generate_figure2(dataset='hg00731', output_suffix='', panel_c_style: str = "summary_bars"):
    """Generate complete Figure 2."""

    setup_style()

    # Use constrained_layout for automatic spacing of labels, titles, and legends
    # Nature Methods maximum width: 180mm = 7.087 inches (use 6.9 to account for bbox_inches='tight')
    fig = plt.figure(figsize=(6.9, 3.1), layout="constrained")
    outer = fig.add_gridspec(1, 3, width_ratios=[0.85, 2.0, 0.95], wspace=0.10)

    ax_a = fig.add_subplot(outer[0, 0])

    inner_b = outer[0, 1].subgridspec(1, 3, wspace=0.25)
    ax_b1 = fig.add_subplot(inner_b[0, 0])
    ax_b2 = fig.add_subplot(inner_b[0, 1])
    ax_b3 = fig.add_subplot(inner_b[0, 2])

    inner_c = outer[0, 2].subgridspec(2, 1, hspace=0.30)
    ax_c1 = fig.add_subplot(inner_c[0, 0])
    ax_c2 = fig.add_subplot(inner_c[1, 0])

    panel_a_speed_comparison(ax_a, dataset=dataset)
    panel_b_count_comparison_3way([ax_b1, ax_b2, ax_b3], dataset=dataset)
    if panel_c_style == "histograms":
        # New default: clear before/after histograms with n>=10 filter (per user feedback)
        panel_c_histograms(ax_c1, ax_c2, dataset=dataset, min_reads=10)
    elif panel_c_style == "scatter_stratified":
        # Scatter plot + stratified bias reduction bars
        panel_c_scatter_stratified(ax_c1, ax_c2, dataset=dataset, min_reads=10)
    elif panel_c_style in ("summary_bars", "summary"):
        panel_c_summary_bars(ax_c1, ax_c2, dataset=dataset)
    elif panel_c_style in ("ref_ratio_hist", "ref_ratio"):
        panel_c_ref_ratio_before_after(ax_c1, ax_c2, dataset=dataset)
    elif panel_c_style in ("two_panel", "bias_two"):
        panel_c_two_plots(ax_c1, ax_c2, dataset=dataset)
    elif panel_c_style == "bias_reduction":
        panel_c_bias_reduction(ax_c1, dataset=dataset)
        ax_c2.axis("off")
    elif panel_c_style == "mean_ci":
        panel_c_delta_mean_ci(ax_c1, dataset=dataset)
        ax_c2.axis("off")
    elif panel_c_style == "delta_skew":
        panel_c_delta_skew(ax_c1, dataset=dataset)
        ax_c2.axis("off")
    elif panel_c_style == "affected_only":
        panel_c_delta_counts_affected_only(ax_c1, dataset=dataset)
        ax_c2.axis("off")
    else:
        panel_c_delta_counts(ax_c1, dataset=dataset)
        ax_c2.axis("off")

    # Title - shorter for constrained_layout
    dataset_name = 'RNA-seq' if dataset == 'hg00731' else 'ATAC-seq'
    fig.suptitle(f'Figure 2: Allele Counting - {dataset_name}',
                 fontsize=11, fontweight='bold')

    # constrained_layout handles spacing automatically - no manual adjustments needed

    # Save
    out_dir = get_plot_path(2, 'figure2').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Add suffix if provided
    base_name = f'figure2{output_suffix}' if output_suffix else 'figure2'

    png_path = out_dir / f'{base_name}.png'
    pdf_path = out_dir / f'{base_name}.pdf'

    plt.savefig(png_path, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    return fig


def generate_combined_figure():
    """Generate combined figure with both datasets."""

    setup_style()

    # Use constrained_layout for automatic spacing
    fig = plt.figure(figsize=(11, 7), layout="constrained")

    # Top row: HG00731
    ax1 = plt.subplot(2, 3, 1)
    ax2 = plt.subplot(2, 3, 2)
    ax3 = plt.subplot(2, 3, 3)

    panel_a_speed_comparison(ax1, dataset='hg00731')
    panel_b_count_comparison(ax2, dataset='hg00731')
    panel_c_bias_reduction(ax3, dataset='hg00731')

    # Bottom row: GM12878
    ax4 = plt.subplot(2, 3, 4)
    ax5 = plt.subplot(2, 3, 5)
    ax6 = plt.subplot(2, 3, 6)

    panel_a_speed_comparison(ax4, dataset='gm12878')
    panel_b_count_comparison(ax5, dataset='gm12878')
    panel_c_bias_reduction(ax6, dataset='gm12878')

    # Add row labels
    ax1.text(-0.25, 0.5, 'RNA-seq\n(HG00731)', transform=ax1.transAxes,
             fontsize=10, fontweight='bold', va='center', ha='right', rotation=90)
    ax4.text(-0.25, 0.5, 'ATAC-seq\n(GM12878)', transform=ax4.transAxes,
             fontsize=10, fontweight='bold', va='center', ha='right', rotation=90)

    # Relabel panels
    for i, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6]):
        label = chr(65 + i)  # A, B, C, D, E, F
        ax.set_title(label, fontsize=9, fontweight='bold', loc='left', x=-0.12)

    fig.suptitle('Figure 2: Allele Counting',
                 fontsize=11, fontweight='bold')

    # constrained_layout handles spacing automatically

    # Save
    out_dir = get_plot_path(2, 'figure2').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    png_path = out_dir / 'figure2_combined.png'
    pdf_path = out_dir / 'figure2_combined.pdf'

    plt.savefig(png_path, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')

    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")

    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Generate Figure 2 plots'
    )
    parser.add_argument(
        '--dataset',
        choices=['hg00731', 'gm12878', 'both', 'combined'],
        default='hg00731',
        help='Dataset to plot (default: hg00731)'
    )
    parser.add_argument(
        '--panel-c-style',
        choices=['histograms', 'scatter_stratified', 'summary_bars', 'ref_ratio_hist', 'two_panel', 'bias_reduction', 'mean_ci', 'delta_skew', 'affected_only', 'all'],
        default='histograms',
        help='Panel C visualization style (default: histograms). "scatter_stratified" shows hexbin scatter + stratified bias reduction bars.'
    )

    args = parser.parse_args()

    if args.panel_c_style == "all":
        # Write all variants for review
        for style, suffix in [
            ("histograms", "_c_histograms"),  # New default (per user feedback)
            ("scatter_stratified", "_c_scatter_stratified"),
            ("summary_bars", "_c_summary_bars"),
            ("ref_ratio_hist", "_c_ref_ratio_hist"),
            ("two_panel", "_c_two_panel"),
            ("bias_reduction", "_c_bias_reduction"),
            ("mean_ci", "_c_mean_ci"),
            ("delta_skew", "_c_delta_skew"),
            ("affected_only", "_c_affected_only"),
        ]:
            if args.dataset in ("hg00731", "gm12878"):
                generate_figure2(dataset=args.dataset, output_suffix=suffix, panel_c_style=style)
            elif args.dataset == "both":
                generate_figure2(dataset='hg00731', output_suffix=f"_hg00731{suffix}", panel_c_style=style)
                generate_figure2(dataset='gm12878', output_suffix=f"_gm12878{suffix}", panel_c_style=style)
            elif args.dataset == "combined":
                # Combined layout not supported for multi-style export; fall back to hg00731.
                generate_figure2(dataset='hg00731', output_suffix=suffix, panel_c_style=style)
        return

    if args.dataset == 'combined':
        # Generate single figure with both datasets
        generate_combined_figure()

    elif args.dataset == 'both':
        # Generate separate figures for each dataset
        print("Generating HG00731 figure...")
        generate_figure2(dataset='hg00731', output_suffix='_hg00731')

        print("\nGenerating GM12878 figure...")
        generate_figure2(dataset='gm12878', output_suffix='_gm12878')

    else:
        # Generate figure for single dataset
        generate_figure2(dataset=args.dataset, panel_c_style=args.panel_c_style)


if __name__ == '__main__':
    main()
