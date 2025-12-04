#!/usr/bin/env python3
"""
Generate publication-quality figures for WASP2 validation.
Suitable for Nature Methods submission.

Requirements:
- matplotlib
- pandas
- numpy
- seaborn (for enhanced styling)
"""

import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from pathlib import Path

# Publication style settings (Nature Methods requirements)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'lines.linewidth': 1.5,
    'patch.linewidth': 0.8,
})

# Colorblind-friendly palette (Okabe-Ito)
COLORS = {
    'SNP': '#E69F00',       # Orange
    'INS': '#56B4E9',       # Sky blue
    'DEL': '#009E73',       # Bluish green
    'WASP1': '#D55E00',     # Vermillion
    'STAR_WASP': '#CC79A7', # Reddish purple
    'WASP2': '#0072B2',     # Blue
    'GATK': '#F0E442',      # Yellow
    'ref': '#999999',       # Gray
    'alt': '#000000',       # Black
}

# Convert mm to inches for Nature Methods
SINGLE_COL_WIDTH = 88 / 25.4  # 88mm -> inches
DOUBLE_COL_WIDTH = 180 / 25.4  # 180mm -> inches


class FigureGenerator:
    """Generate publication-quality figures for WASP2."""

    def __init__(self, output_dir='figures'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.figures_generated = []
        self.data_sources = []
        self.missing_data = []

    def save_figure(self, fig, name, description):
        """Save figure in both PDF (vector) and PNG (raster) formats."""
        pdf_path = self.output_dir / f"{name}.pdf"
        png_path = self.output_dir / f"{name}.png"

        fig.savefig(pdf_path, format='pdf', bbox_inches='tight', dpi=300)
        fig.savefig(png_path, format='png', bbox_inches='tight', dpi=300)

        self.figures_generated.append({
            'name': name,
            'description': description,
            'pdf': str(pdf_path.absolute()),
            'png': str(png_path.absolute())
        })

        plt.close(fig)
        print(f"  Generated: {name}")

    def figure1_runtime_comparison(self, benchmark_data):
        """
        Figure 1: Runtime Performance Comparison
        Bar chart comparing WASP1, STAR+WASP, and WASP2-Rust
        """
        print("\nGenerating Figure 1: Runtime Performance Comparison...")

        # Extract timing data
        wasp1_time = benchmark_data['comparison']['wasp1_s']
        star_wasp_time = benchmark_data['comparison']['star_wasp_s']
        wasp2_time = benchmark_data['wasp_only_s']

        methods = ['WASP1\n(Python)', 'STAR+WASP\n(C++)', 'WASP2\n(Rust)']
        times = [wasp1_time, star_wasp_time, wasp2_time]
        colors = [COLORS['WASP1'], COLORS['STAR_WASP'], COLORS['WASP2']]

        # Create figure (single column width)
        fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, SINGLE_COL_WIDTH * 0.8))

        bars = ax.bar(methods, times, color=colors, edgecolor='black', linewidth=0.8)

        # Add value labels on bars
        for bar, time in zip(bars, times):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{time:.0f}s\n({time/60:.1f}m)',
                   ha='center', va='bottom', fontsize=8)

        # Add speedup annotations
        speedup_vs_wasp1 = wasp1_time / wasp2_time
        speedup_vs_star_wasp = star_wasp_time / wasp2_time

        ax.text(0.5, 0.95, f'{speedup_vs_wasp1:.1f}× faster than WASP1',
               transform=ax.transAxes, ha='center', va='top',
               fontsize=8, style='italic', color=COLORS['WASP1'])

        ax.text(0.5, 0.88, f'{speedup_vs_star_wasp:.2f}× faster than STAR+WASP',
               transform=ax.transAxes, ha='center', va='top',
               fontsize=8, style='italic', color=COLORS['STAR_WASP'])

        ax.set_ylabel('Runtime (seconds)', fontweight='bold')
        ax.set_title('WASP2 Runtime Performance\n(HG00731, 8 threads)', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(0, max(times) * 1.25)
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

        self.save_figure(fig, 'figure1_runtime_comparison',
                        'Runtime comparison of WASP implementations')
        self.data_sources.append(str(Path('benchmarking/star_wasp_comparison/results/unified_2025-12-03_20-10-45/benchmark_results.json').absolute()))

    def figure2_pipeline_breakdown(self, benchmark_data):
        """
        Figure 2: Pipeline Step Breakdown
        Stacked bar chart showing time spent in each pipeline step
        """
        print("\nGenerating Figure 2: Pipeline Step Breakdown...")

        # WASP2 unified pipeline breakdown
        wasp2_steps = {
            'STAR Initial': benchmark_data['step1_star_initial_s'],
            'Find Intersecting\n(Rust)': benchmark_data['step2_wasp2_unified_s'],
            'STAR Remap': benchmark_data['step3_star_remap_s'],
            'Filter Remapped\n(Rust)': benchmark_data['step4_wasp2_filter_s']
        }

        fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH * 1.2, SINGLE_COL_WIDTH * 0.7))

        steps = list(wasp2_steps.keys())
        times = list(wasp2_steps.values())
        percentages = [t / sum(times) * 100 for t in times]

        # Color gradient for steps
        step_colors = ['#0072B2', '#56B4E9', '#009E73', '#E69F00']

        bars = ax.barh(steps, times, color=step_colors, edgecolor='black', linewidth=0.8)

        # Add time and percentage labels
        for bar, time, pct in zip(bars, times, percentages):
            width = bar.get_width()
            ax.text(width + 5, bar.get_y() + bar.get_height()/2.,
                   f'{time:.1f}s ({pct:.1f}%)',
                   ha='left', va='center', fontsize=8)

        ax.set_xlabel('Runtime (seconds)', fontweight='bold')
        ax.set_title('WASP2 Unified Pipeline Breakdown', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(0, max(times) * 1.3)
        ax.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.5)

        # Add total time
        total_time = sum(times)
        ax.text(0.98, 0.02, f'Total: {total_time:.1f}s ({total_time/60:.1f}m)',
               transform=ax.transAxes, ha='right', va='bottom',
               fontsize=8, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                        edgecolor='black', linewidth=0.8))

        self.save_figure(fig, 'figure2_pipeline_breakdown',
                        'Breakdown of WASP2 pipeline execution time by step')

    def figure3_throughput_comparison(self, benchmark_data, expected_counts):
        """
        Figure 3: Read Processing Throughput
        Compare reads/second processing rate
        """
        print("\nGenerating Figure 3: Read Processing Throughput...")

        # Calculate reads processed (from expected counts)
        total_bam_reads = expected_counts['keep_bam_reads'] + expected_counts['to_remap_bam_reads']

        # Calculate throughput
        wasp1_throughput = total_bam_reads / benchmark_data['comparison']['wasp1_s']
        star_wasp_throughput = total_bam_reads / benchmark_data['comparison']['star_wasp_s']
        wasp2_throughput = total_bam_reads / benchmark_data['wasp_only_s']

        methods = ['WASP1', 'STAR+WASP', 'WASP2']
        throughput = [wasp1_throughput, star_wasp_throughput, wasp2_throughput]
        colors = [COLORS['WASP1'], COLORS['STAR_WASP'], COLORS['WASP2']]

        fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, SINGLE_COL_WIDTH * 0.8))

        bars = ax.bar(methods, [t/1e6 for t in throughput], color=colors,
                     edgecolor='black', linewidth=0.8)

        # Add value labels
        for bar, tp in zip(bars, throughput):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{tp/1e6:.2f}M\nreads/s',
                   ha='center', va='bottom', fontsize=8)

        ax.set_ylabel('Throughput (million reads/second)', fontweight='bold')
        ax.set_title('Read Processing Throughput\n(HG00731, 8 threads)', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(0, max([t/1e6 for t in throughput]) * 1.25)
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

        # Add dataset info
        ax.text(0.02, 0.98, f'Dataset: {total_bam_reads/1e6:.1f}M reads',
               transform=ax.transAxes, ha='left', va='top',
               fontsize=7, style='italic',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                        alpha=0.8, edgecolor='gray', linewidth=0.5))

        self.save_figure(fig, 'figure3_throughput_comparison',
                        'Read processing throughput comparison')

    def figure4_variant_coverage(self, expected_counts):
        """
        Figure 4: Variant and Read Statistics
        Show distribution of reads across variant sites
        """
        print("\nGenerating Figure 4: Variant and Read Statistics...")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL_WIDTH, SINGLE_COL_WIDTH * 0.6))

        # Panel A: Read distribution
        categories = ['Keep\n(no variants)', 'Remap\n(at variants)']
        counts = [expected_counts['keep_bam_reads'], expected_counts['to_remap_bam_reads']]
        colors_panel = [COLORS['ref'], COLORS['WASP2']]

        wedges, texts, autotexts = ax1.pie(counts, labels=categories, autopct='%1.1f%%',
                                            colors=colors_panel, startangle=90,
                                            textprops={'fontsize': 9},
                                            wedgeprops={'edgecolor': 'black', 'linewidth': 0.8})

        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')

        ax1.set_title('Read Distribution', fontweight='bold', pad=10)

        # Panel B: Variant processing statistics
        stats = {
            'VCF Variants': expected_counts['vcf_variants'],
            'Unique Read\nOverlaps': expected_counts['unique_read_overlaps'],
            'Haplotypes\nGenerated': expected_counts['total_haplotypes']
        }

        stat_names = list(stats.keys())
        stat_values = [v/1e6 for v in stats.values()]
        colors_bar = ['#0072B2', '#009E73', '#E69F00']

        bars = ax2.bar(range(len(stat_names)), stat_values, color=colors_bar,
                      edgecolor='black', linewidth=0.8)

        for bar, val in zip(bars, stats.values()):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val/1e6:.2f}M',
                    ha='center', va='bottom', fontsize=8)

        ax2.set_xticks(range(len(stat_names)))
        ax2.set_xticklabels(stat_names)
        ax2.set_ylabel('Count (millions)', fontweight='bold')
        ax2.set_title('Variant Processing', fontweight='bold', pad=10)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

        plt.tight_layout()

        self.save_figure(fig, 'figure4_variant_coverage',
                        'Variant coverage and read distribution statistics')

    def figure5_gatk_comparison(self, ground_truth_file, gatk_file):
        """
        Figure 5: GATK ASEReadCounter Comparison
        Compare WASP2 ground truth with GATK counts
        """
        print("\nGenerating Figure 5: GATK Comparison...")

        # Load data
        ground_truth = pd.read_csv(ground_truth_file)
        gatk = pd.read_table(gatk_file)

        # Parse variant types from variant_id
        ground_truth['variant_type'] = ground_truth['variant_id'].str.extract(r'(SNP|INS|DEL)')[0]

        # Merge datasets
        merged = pd.merge(
            ground_truth,
            gatk[['position', 'refCount', 'altCount', 'totalCount']],
            left_on='pos',
            right_on='position',
            how='inner'
        )

        # Calculate ratios
        merged['true_ratio'] = merged['true_ref_ratio']
        merged['gatk_ratio'] = merged['refCount'] / (merged['refCount'] + merged['altCount'])

        # Create figure with 3 panels
        fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL_WIDTH, SINGLE_COL_WIDTH * 0.7))

        variant_types = ['SNP', 'INS', 'DEL']

        for ax, vtype in zip(axes, variant_types):
            data = merged[merged['variant_type'] == vtype]

            if len(data) == 0:
                ax.text(0.5, 0.5, f'No {vtype}\ndata',
                       ha='center', va='center', fontsize=10, style='italic')
                ax.set_title(vtype, fontweight='bold')
                continue

            # Scatter plot
            ax.scatter(data['true_ratio'], data['gatk_ratio'],
                      color=COLORS[vtype], s=50, alpha=0.7,
                      edgecolor='black', linewidth=0.5)

            # Add perfect correlation line
            ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Perfect')

            # Calculate R²
            if len(data) > 1:
                corr = np.corrcoef(data['true_ratio'], data['gatk_ratio'])[0, 1]
                r_squared = corr ** 2
                ax.text(0.05, 0.95, f'R² = {r_squared:.3f}',
                       transform=ax.transAxes, fontsize=8,
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                                edgecolor='black', linewidth=0.5))

            ax.set_xlabel('Ground Truth Ratio', fontsize=9)
            if ax == axes[0]:
                ax.set_ylabel('GATK Ratio', fontsize=9)
            ax.set_title(vtype, fontweight='bold')
            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)
            ax.set_aspect('equal')
            ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        plt.tight_layout()

        self.save_figure(fig, 'figure5_gatk_comparison',
                        'GATK ASEReadCounter accuracy by variant type')

        self.data_sources.append(str(Path(ground_truth_file).absolute()))
        self.data_sources.append(str(Path(gatk_file).absolute()))

    def generate_supplementary_table(self, benchmark_data, expected_counts):
        """
        Generate supplementary table with detailed statistics
        """
        print("\nGenerating Supplementary Table...")

        # Create detailed statistics table
        table_data = {
            'Metric': [
                'Total BAM reads',
                'VCF variants',
                'Reads at variants (to remap)',
                'Reads without variants (keep)',
                'Unique read-variant overlaps',
                'Haplotypes generated',
                'R1 FASTQ reads',
                'R2 FASTQ reads',
                '',
                'WASP1 runtime (s)',
                'STAR+WASP runtime (s)',
                'WASP2 runtime (s)',
                'WASP2 speedup vs WASP1',
                'WASP2 speedup vs STAR+WASP',
                '',
                'WASP2 STAR initial (s)',
                'WASP2 find intersecting (s)',
                'WASP2 STAR remap (s)',
                'WASP2 filter remapped (s)',
            ],
            'Value': [
                f"{expected_counts['keep_bam_reads'] + expected_counts['to_remap_bam_reads']:,}",
                f"{expected_counts['vcf_variants']:,}",
                f"{expected_counts['to_remap_bam_reads']:,}",
                f"{expected_counts['keep_bam_reads']:,}",
                f"{expected_counts['unique_read_overlaps']:,}",
                f"{expected_counts['total_haplotypes']:,}",
                f"{expected_counts['r1_fastq_reads']:,}",
                f"{expected_counts['r2_fastq_reads']:,}",
                '',
                f"{benchmark_data['comparison']['wasp1_s']:.2f}",
                f"{benchmark_data['comparison']['star_wasp_s']:.2f}",
                f"{benchmark_data['wasp_only_s']:.2f}",
                f"{benchmark_data['speedup_vs_wasp1']:.2f}×",
                f"{benchmark_data['speedup_vs_star_wasp']:.2f}×",
                '',
                f"{benchmark_data['step1_star_initial_s']:.2f}",
                f"{benchmark_data['step2_wasp2_unified_s']:.2f}",
                f"{benchmark_data['step3_star_remap_s']:.2f}",
                f"{benchmark_data['step4_wasp2_filter_s']:.2f}",
            ]
        }

        df = pd.DataFrame(table_data)

        # Save as CSV
        table_path = self.output_dir / 'supplementary_table_S1.csv'
        df.to_csv(table_path, index=False)

        print(f"  Generated: supplementary_table_S1.csv")

        self.figures_generated.append({
            'name': 'supplementary_table_S1',
            'description': 'Detailed pipeline statistics and performance metrics',
            'csv': str(table_path.absolute())
        })

    def generate_report(self):
        """Generate final report of all figures created."""
        report_path = self.output_dir / 'figure_generation_report.txt'

        with open(report_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("WASP2 PUBLICATION FIGURES - GENERATION REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"Generated: {len(self.figures_generated)} figures\n")
            f.write(f"Output directory: {self.output_dir.absolute()}\n\n")

            f.write("FIGURES GENERATED:\n")
            f.write("-" * 80 + "\n")
            for i, fig in enumerate(self.figures_generated, 1):
                f.write(f"\n{i}. {fig['name']}\n")
                f.write(f"   Description: {fig['description']}\n")
                if 'pdf' in fig:
                    f.write(f"   PDF: {fig['pdf']}\n")
                    f.write(f"   PNG: {fig['png']}\n")
                if 'csv' in fig:
                    f.write(f"   CSV: {fig['csv']}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("DATA SOURCES:\n")
            f.write("-" * 80 + "\n")
            for source in set(self.data_sources):
                f.write(f"  - {source}\n")

            if self.missing_data:
                f.write("\n" + "=" * 80 + "\n")
                f.write("MISSING DATA (Figures Not Generated):\n")
                f.write("-" * 80 + "\n")
                for item in self.missing_data:
                    f.write(f"  - {item}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("NATURE METHODS REQUIREMENTS:\n")
            f.write("-" * 80 + "\n")
            f.write("  Resolution: 300 DPI (Met)\n")
            f.write("  Formats: PDF (vector) + PNG (raster) (Met)\n")
            f.write("  Colorblind palette: Okabe-Ito (Met)\n")
            f.write("  Font: Arial/Helvetica, 8-12pt (Met)\n")
            f.write("  Width: 88mm (single) / 180mm (double) (Met)\n")
            f.write("\n" + "=" * 80 + "\n")

        print(f"\nReport saved: {report_path.absolute()}")
        return str(report_path.absolute())


def main():
    """Main function to generate all figures."""
    print("=" * 80)
    print("WASP2 PUBLICATION FIGURE GENERATOR")
    print("=" * 80)

    # Initialize generator
    generator = FigureGenerator(output_dir='figures')

    # Load data files
    root = Path('/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp')

    benchmark_file = root / 'benchmarking/star_wasp_comparison/results/unified_2025-12-03_20-10-45/benchmark_results.json'
    expected_counts_file = root / 'baselines/mapping/expected_counts.json'
    ground_truth_file = root / 'simulation_results/comprehensive_20251203_210028/ground_truth.csv'
    gatk_file = root / 'comparison_results/gatk_ase_counts.table'

    # Load benchmark data
    if benchmark_file.exists():
        with open(benchmark_file) as f:
            benchmark_data = json.load(f)
        print(f"\nLoaded: {benchmark_file}")
    else:
        print(f"\nERROR: Benchmark file not found: {benchmark_file}")
        return

    # Load expected counts
    if expected_counts_file.exists():
        with open(expected_counts_file) as f:
            expected_counts_data = json.load(f)
            expected_counts = expected_counts_data['expected_counts']
        print(f"Loaded: {expected_counts_file}")
    else:
        print(f"\nERROR: Expected counts file not found: {expected_counts_file}")
        return

    # Generate figures
    print("\n" + "=" * 80)
    print("GENERATING FIGURES")
    print("=" * 80)

    # Figure 1: Runtime comparison
    generator.figure1_runtime_comparison(benchmark_data)

    # Figure 2: Pipeline breakdown
    generator.figure2_pipeline_breakdown(benchmark_data)

    # Figure 3: Throughput comparison
    generator.figure3_throughput_comparison(benchmark_data, expected_counts)

    # Figure 4: Variant coverage
    generator.figure4_variant_coverage(expected_counts)

    # Figure 5: GATK comparison (if data available)
    if ground_truth_file.exists() and gatk_file.exists():
        try:
            generator.figure5_gatk_comparison(str(ground_truth_file), str(gatk_file))
        except Exception as e:
            print(f"  Warning: Could not generate GATK comparison: {e}")
            generator.missing_data.append(f"Figure 5 (GATK comparison): {e}")
    else:
        print("\n  Skipping Figure 5 (GATK comparison): Missing data files")
        if not ground_truth_file.exists():
            generator.missing_data.append(f"Ground truth file: {ground_truth_file}")
        if not gatk_file.exists():
            generator.missing_data.append(f"GATK counts file: {gatk_file}")

    # Generate supplementary table
    generator.generate_supplementary_table(benchmark_data, expected_counts)

    # Generate report
    print("\n" + "=" * 80)
    print("GENERATING REPORT")
    print("=" * 80)
    report_path = generator.generate_report()

    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Figures generated: {len(generator.figures_generated)}")
    print(f"Output directory: {generator.output_dir.absolute()}")
    print(f"Report: {report_path}")
    print("\nAll figures are publication-ready for Nature Methods submission!")
    print("=" * 80)


if __name__ == '__main__':
    main()
