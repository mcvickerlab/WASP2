# WASP2 Publication Figures

Publication-quality figures for WASP2 validation, prepared for **Nature Methods** submission.

## Generation Information

- **Generated:** December 3, 2025
- **Script:** `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/scripts/generate_publication_figures.py`
- **Total Figures:** 5 main figures + 1 supplementary table

## Nature Methods Compliance

All figures meet Nature Methods requirements:

- **Resolution:** 300 DPI (both PDF and PNG)
- **Formats:** PDF (vector graphics) + PNG (raster backup)
- **Color Palette:** Colorblind-friendly (Okabe-Ito palette)
- **Typography:** Arial/Helvetica, 8-12pt
- **Dimensions:**
  - Single column: 88mm (3.46")
  - Double column: 180mm (7.09")

## Figures

### Figure 1: Runtime Performance Comparison
**File:** `figure1_runtime_comparison.pdf/png`
**Description:** Bar chart comparing runtime performance of WASP1 (Python), STAR+WASP (C++), and WASP2 (Rust)
**Key Results:**
- WASP2: 29.85s
- STAR+WASP: 331.67s (0.59× slower than WASP2)
- WASP1: 1541.13s (2.77× slower than WASP2)
- **Layout:** Single column width (88mm)

### Figure 2: Pipeline Step Breakdown
**File:** `figure2_pipeline_breakdown.pdf/png`
**Description:** Horizontal bar chart showing time distribution across WASP2 pipeline steps
**Pipeline Steps:**
1. STAR Initial Alignment: 475.95s (85.6%)
2. Find Intersecting Variants (Rust): 16.88s (3.0%)
3. STAR Remap: 48.03s (8.6%)
4. Filter Remapped (Rust): 12.97s (2.3%)
- **Total:** 553.83s
- **Layout:** Single column width (88mm)

### Figure 3: Read Processing Throughput
**File:** `figure3_throughput_comparison.pdf/png`
**Description:** Bar chart comparing read processing throughput (million reads/second)
**Throughput Results:**
- WASP2: 3.77M reads/s
- STAR+WASP: 0.34M reads/s
- WASP1: 0.07M reads/s
- **Dataset:** 112.5M reads from HG00731
- **Layout:** Single column width (88mm)

### Figure 4: Variant and Read Statistics
**File:** `figure4_variant_coverage.pdf/png`
**Description:** Two-panel figure showing read distribution and variant processing statistics
**Panel A (Pie Chart):** Read Distribution
- Reads without variants (keep): 99.9M (88.8%)
- Reads at variants (remap): 12.6M (11.2%)

**Panel B (Bar Chart):** Variant Processing
- VCF Variants: 2.21M
- Unique Read-Variant Overlaps: 6.89M
- Haplotypes Generated: 3.07M
- **Layout:** Double column width (180mm)

### Figure 5: GATK ASEReadCounter Comparison
**File:** `figure5_gatk_comparison.pdf/png`
**Description:** Three-panel scatter plot comparing ground truth allelic ratios with GATK ASEReadCounter predictions
**Panels:**
- SNP: Correlation with ground truth
- INS (Insertion): Correlation with ground truth
- DEL (Deletion): Correlation with ground truth

**Purpose:** Demonstrates GATK's handling of different variant types compared to known ground truth from simulation
- **Layout:** Double column width (180mm)

## Supplementary Materials

### Supplementary Table S1
**File:** `supplementary_table_S1.csv`
**Description:** Detailed pipeline statistics and performance metrics
**Contents:**
- Read counts (total, remapped, kept)
- Variant statistics
- Runtime comparisons (WASP1, STAR+WASP, WASP2)
- Pipeline step timings

**Format:** CSV (comma-separated values) for easy import into LaTeX tables

## Data Sources

All figures were generated from the following validated data sources:

1. **Benchmark Data:**
   - `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/star_wasp_comparison/results/unified_2025-12-03_20-10-45/benchmark_results.json`

2. **Expected Counts (Validation Baseline):**
   - `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/baselines/mapping/expected_counts.json`

3. **Simulation Ground Truth:**
   - `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/simulation_results/comprehensive_20251203_210028/ground_truth.csv`

4. **GATK Comparison:**
   - `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/comparison_results/gatk_ase_counts.table`

## Dataset Details

**Sample:** HG00731 (1000 Genomes Project)
**Data Type:** RNA-seq, paired-end
**Total Reads:** 112,509,428
**VCF Variants:** 2,206,785 heterozygous sites
**Processing Environment:**
- CPU: 8 threads
- Platform: Linux x86_64
- Branch: `ropc-indels` (INDEL support)

## Color Palette (Okabe-Ito)

Colorblind-friendly colors used throughout:

- **SNP:** #E69F00 (Orange)
- **INS:** #56B4E9 (Sky Blue)
- **DEL:** #009E73 (Bluish Green)
- **WASP1:** #D55E00 (Vermillion)
- **STAR+WASP:** #CC79A7 (Reddish Purple)
- **WASP2:** #0072B2 (Blue)
- **GATK:** #F0E442 (Yellow)

## Reproduction

To regenerate these figures:

```bash
# Activate environment
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2

# Navigate to project directory
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

# Run figure generation script
python scripts/generate_publication_figures.py
```

## Dependencies

Python packages required:
- matplotlib >= 3.5.0
- pandas >= 1.4.0
- numpy >= 1.21.0

All dependencies are installed in the `WASP2_dev2` conda environment.

## File Inventory

```
figures/
├── figure1_runtime_comparison.pdf        (31 KB, vector)
├── figure1_runtime_comparison.png        (112 KB, 300 DPI)
├── figure2_pipeline_breakdown.pdf        (27 KB, vector)
├── figure2_pipeline_breakdown.png        (101 KB, 300 DPI)
├── figure3_throughput_comparison.pdf     (27 KB, vector)
├── figure3_throughput_comparison.png     (83 KB, 300 DPI)
├── figure4_variant_coverage.pdf          (26 KB, vector)
├── figure4_variant_coverage.png          (104 KB, 300 DPI)
├── figure5_gatk_comparison.pdf           (17 KB, vector)
├── figure5_gatk_comparison.png           (80 KB, 300 DPI)
├── supplementary_table_S1.csv            (551 B)
└── figure_generation_report.txt          (3.8 KB)
```

## Notes

- PDF files are recommended for publication as they contain vector graphics
- PNG files are provided as high-resolution raster backups (300 DPI)
- All figures use consistent styling and color schemes
- Figures are ready for direct inclusion in manuscript LaTeX source

## Citation

When using these figures, please cite:

> WASP2: Rust-Accelerated Allele-Specific Read Mapping
> [Authors]
> Nature Methods (in preparation)

---

**Last Updated:** December 3, 2025
**Contact:** jjaureguy@anthropic.com
