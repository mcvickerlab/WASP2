"""
WASP2 Paper - Shared Configuration

All paths and settings for figure generation.
"""
from pathlib import Path

# =============================================================================
# BASE PATHS
# =============================================================================
REPO_ROOT = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp")
PAPER_DIR = REPO_ROOT / "paper"

# =============================================================================
# INPUT DATA PATHS
# =============================================================================

# Reference genomes
REF_GENOME_HG38 = Path("/iblm/netapp/data1/aho/ref_genomes/index/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa")
STAR_INDEX = Path("/iblm/netapp/data1/external/GRC38/combined/GDC_hg38/star_index")

# Variant files
VCF_NA12878 = Path("/iblm/netapp/data1/aho/variants/NA12878.vcf.gz")
VCF_HG00731 = REPO_ROOT / "benchmarking/star_wasp_comparison/data/HG00731_het_only_chr.vcf.gz"

# BAM files - GM12878 ATAC-seq
BAM_GM12878_ATAC_RAW = Path("/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam")

# BAM files - HG00731 RNA-seq
BAM_HG00731_RNA = REPO_ROOT / "benchmarking/star_wasp_comparison/results/unified_2025-12-04_00-29-39/A_sorted.bam"

# Gene annotations
IMPRINTED_GENES = Path("/iblm/netapp/home/aho/projects/wasp/testing/performance/data/GM12878_ATACseq_50k_merged/test_imprinted/geneimprint_coords.tsv")

# =============================================================================
# EXISTING RESULT FILES
# =============================================================================

# Gene imprinting results
GENE_IMPRINTING_RESULTS_DIR = REPO_ROOT / "validation/gene_imprinting/results/rust_pipeline"
GENE_SIGNIFICANCE_TSV = GENE_IMPRINTING_RESULTS_DIR / "gene_significance.tsv"
GENE_IMPRINTING_FIGURE = GENE_IMPRINTING_RESULTS_DIR / "figure_publication.png"

# Benchmark results
BENCHMARK_RESULTS_DIR = REPO_ROOT / "benchmarking/results"
BENCHMARK_4PANEL_FIGURE = REPO_ROOT / "plots/comprehensive_benchmark_4panel.png"
PHASER_1THREAD_TIMING = BENCHMARK_RESULTS_DIR / "phaser_1thread/timing.csv"

# =============================================================================
# TOOL PATHS
# =============================================================================
BWA = Path("/iblm/netapp/oldhome/j3gu/anaconda2/bin/bwa")
SAMTOOLS = Path("/iblm/netapp/home/jjaureguy/mambaforge/bin/samtools")
GATK = Path("/iblm/netapp/home/jjaureguy/mambaforge/envs/WASP2_dev2/bin/gatk")

# =============================================================================
# BENCHMARK DATA (from previous runs)
# =============================================================================
BENCHMARK_DATA = {
    # Thread scaling - WASP2-Rust (allele counting)
    'wasp2_rust_scaling': {
        1: {'time': 57.99, 'variants_per_sec': 38058},
        2: {'time': 36.28, 'variants_per_sec': 60829},
        4: {'time': 17.40, 'variants_per_sec': 126844},
        8: {'time': 9.32, 'variants_per_sec': 236677},
    },
    # phASER (allele counting)
    'phaser': {
        1: {'time': 561.0},
        8: {'time': 162.0},
    },
    # GATK ASEReadCounter (allele counting)
    'gatk': {
        1: {'time': 1600.0},
    },
    # WASP versions comparison (at 150M reads, 8 threads)
    'wasp_versions': {
        'wasp1': {'time': 7000},
        'wasp2_python': {'time': 2000},
        'wasp2_rust_snp': {'time': 700},
        'wasp2_rust_indel': {'time': 750},
    },
    # ==========================================================================
    # ATAC-seq Read Filtering Benchmark (GM12878, 159M reads, 8 threads)
    # ==========================================================================
    'atac_filtering': {
        'sample': 'NA12878',
        'assay': 'ATAC-seq',
        'reads': 159_133_683,
        'threads': 8,
        'aligner': 'BWA-MEM',
        'wasp1': {
            'wasp_only_s': 4816.67,
            'total_s': 5523.97,
        },
        'wasp2_rust_snp': {
            'wasp_only_s': 78.88,
            'total_s': 996.40,
        },
        'wasp2_rust_indel': {
            'wasp_only_s': 77.95,
            'total_s': 626.20,
        },
    },
    # ==========================================================================
    # RNA-seq Read Filtering Benchmark (HG00731, 56M reads, 8 threads)
    # Source: benchmarking/star_wasp_comparison/results/
    # ==========================================================================
    'rnaseq_filtering': {
        'sample': 'HG00731',
        'assay': 'RNA-seq',
        'reads': 56_200_000,  # 56.2M paired-end reads (from STAR metrics)
        'threads': 8,
        'aligner': 'STAR',
        # WASP1: Our actual benchmark run (Dec 6, 2025)
        # File: wasp1_2025-12-06_20-28-15/benchmark_results.json
        'wasp1': {
            'wasp_only_s': 2405.00,  # Steps 2+4+5 (find_intersecting + filter + merge)
            'total_s': 2961.17,  # Full pipeline including STAR alignment
            'source': 'wasp1_2025-12-06_20-28-15/benchmark_results.json',
        },
        # STAR-WASP C++: ACTUAL benchmark run (not published value!)
        # File: star_wasp_2025-12-04_00-58-54/benchmark_results.json
        'star_wasp': {
            'wasp_only_s': None,  # Integrated in aligner
            'total_s': 532.61,  # Our actual benchmark, NOT 331.67s from paper
            'source': 'star_wasp_2025-12-04_00-58-54/benchmark_results.json',
        },
        # WASP2-Rust unified: optimized run
        # File: unified_opt_2025-12-04_02-56-39/benchmark_results.json
        'wasp2_rust': {
            'wasp_only_s': 30.28,
            'total_s': 576.28,
            'source': 'unified_opt_2025-12-04_02-56-39/benchmark_results.json',
        },
    },
}

# Gene imprinting results
IMPRINTING_RESULTS = {
    'snp_only_significant': 2,
    'snp_indel_significant': 5,
    'improvement_factor': 2.5,
    'newly_significant_genes': ['PAOX', 'SNRPN', 'RASGRF1'],
}

# =============================================================================
# OUTPUT PATHS
# =============================================================================
def get_figure_dir(fig_num: int) -> Path:
    """Get directory for a specific figure."""
    return PAPER_DIR / f"figure{fig_num}"

def get_plot_path(fig_num: int, name: str, ext: str = 'png') -> Path:
    """Get path for a plot file."""
    return get_figure_dir(fig_num) / "plots" / f"{name}.{ext}"

def get_data_path(fig_num: int, name: str) -> Path:
    """Get path for a data file."""
    return get_figure_dir(fig_num) / "data" / name

# =============================================================================
# PLOT SETTINGS (Nature Methods style - strict compliance)
# Reference: https://www.nature.com/nmeth/submission-guidelines/aip-and-formatting
# =============================================================================
PLOT_SETTINGS = {
    'font_family': 'sans-serif',
    'font_sans_serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font_size': 7,              # Max 7pt for general text
    'axes_labelsize': 7,         # Axis labels: 7pt
    'axes_titlesize': 8,         # Panel labels: 8pt bold
    'legend_fontsize': 6,        # Legend: 6pt (min 5pt)
    'tick_labelsize': 6,         # Tick labels: 6pt
    'figure_dpi': 150,
    'savefig_dpi': 300,
    # Nature Methods column widths
    'single_column_mm': 88,      # 3.46 inches
    'double_column_mm': 180,     # 7.09 inches
    'one_half_column_mm': 120,   # 4.72 inches
}

# Colorblind-safe palette (Bang Wong, Nature Methods 2011)
COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'yellow': '#F0E442',
    'sky_blue': '#56B4E9',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'black': '#000000',
}

# Tool colors
TOOL_COLORS = {
    'wasp1': COLORS['black'],
    'wasp2_python': COLORS['blue'],
    'wasp2_rust': COLORS['green'],
    'wasp2_rust_indel': COLORS['green'],
    'phaser': COLORS['purple'],
    'gatk': COLORS['vermillion'],
    'star_wasp': COLORS['sky_blue'],
}
