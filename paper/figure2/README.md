# Figure 2: Allele Counting Benchmarks

This directory contains scripts and data for generating Figure 2, which compares WASP2-Rust allele counting performance against GATK ASEReadCounter and phASER.

## Figure 2 Panels

### Panel A: Speed Comparison
Bar chart comparing runtime of WASP2-Rust, GATK ASEReadCounter, and phASER for allele counting on the same dataset. Shows speedup factors.

### Panel B: Count Correlation
3-way count correlation (WASP2 vs GATK, WASP2 vs phASER, GATK vs phASER) on the same processed BAM.

### Panel C: Bias Reduction
Per-site change in **ref/alt counts** after WASP2 processing (Δcounts = processed − original) for each counting method (WASP2-Rust, GATK ASEReadCounter, phASER).

## Datasets

1. **HG00731 RNA-seq**: 56M paired-end reads, STAR-aligned, ~2.2M heterozygous SNPs (primary figure)
2. **GM12878 ATAC-seq**: Optional secondary dataset (supported by scripts, not required for the main RNA-seq figure)

## Directory Structure

```
paper/figure2/
├── README.md                    # This file
├── scripts/
│   ├── run_figure2_benchmarks.sh         # SGE job to run all benchmarks
│   ├── generate_count_comparison.py      # Compare counts across tools
│   ├── generate_bias_comparison.py       # Analyze bias reduction
│   ├── generate_before_after_counts.py   # Build Panel C before/after table (all 3 methods)
│   └── generate_figure2.py               # Generate final plots
├── data/
│   ├── timing_results.json               # Benchmark timing data
│   ├── hg00731/
│   │   ├── wasp2_counts.filtered.tsv     # WASP2-Rust counts (processed BAM)
│   │   ├── wasp2_counts.original.tsv     # WASP2-Rust counts (original BAM)
│   │   ├── gatk_counts.filtered.table    # GATK counts (WASP2-processed BAM)
│   │   ├── gatk_counts.original.table    # GATK counts (original BAM, pre-WASP)
│   │   ├── phaser.filtered.allelic_counts.txt  # phASER counts (processed BAM)
│   │   ├── phaser.original.allelic_counts.txt  # phASER counts (original BAM)
│   │   ├── count_comparison.tsv          # Merged counts
│   │   └── before_after_counts.tsv       # Panel C (before/after across methods)
│   └── gm12878/
│       ├── wasp2_counts.tsv
│       ├── gatk_counts.table
│       ├── phaser.haplotypic_counts.txt
│       ├── count_comparison.tsv
│       └── original_vs_remapped.tsv
├── logs/                                  # Benchmark logs
└── plots/
    ├── figure2.png                        # Main figure (HG00731)
    ├── figure2.pdf
    ├── figure2_hg00731.png                # HG00731-specific
    ├── figure2_gm12878.png                # GM12878-specific
    └── figure2_combined.png               # Both datasets in 2x3 layout
```

## Usage

### Step 1: Run Benchmarks

Submit the SGE job to run all benchmarks:

```bash
cd paper/figure2
mkdir -p logs
qsub scripts/run_figure2_benchmarks.sh
```

This will:
- Run WASP2-Rust allele counting on the selected dataset(s)
- Run GATK ASEReadCounter on both datasets
- Run phASER on both datasets
- Save timing results to `data/timing_results.json`
- Save count outputs to respective dataset directories

**Runtime**: ~2-3 hours on 8 cores, 16GB RAM

**Note**: Make sure the conda environment `WASP2_dev2` is set up with all required tools (GATK, phASER dependencies).

### Step 2: Generate Count Comparisons

After benchmarks complete, compare counts across tools:

```bash
# Process both datasets
python scripts/generate_count_comparison.py --dataset both

# Or process individually
python scripts/generate_count_comparison.py --dataset hg00731
python scripts/generate_count_comparison.py --dataset gm12878
```

This will:
- Parse output from all three tools
- Merge counts by variant position
- Calculate Pearson correlations
- Save merged data to `data/{dataset}/count_comparison.tsv`

### Step 3: Generate Bias Comparisons

Build the Panel C before/after table (all three methods):

```bash
python scripts/generate_before_after_counts.py --dataset hg00731 --min-total 10
python scripts/generate_delta_counts.py --dataset hg00731 --min-total 10
```

This will:
- Merge original vs processed counts for **WASP2, GATK, and phASER**
- Apply a coverage cutoff (`min-total`) consistently
- Save results to `data/hg00731/before_after_counts.tsv`

### Step 4: Generate Plots

Create the final figure:

```bash
# Generate figure for HG00731 only (default)
python scripts/generate_figure2.py

# Generate figure for GM12878 only
python scripts/generate_figure2.py --dataset gm12878

# Generate separate figures for both datasets
python scripts/generate_figure2.py --dataset both

# Generate combined 2x3 figure with both datasets
python scripts/generate_figure2.py --dataset combined
```

Output files will be saved to `plots/`.

## Input Data Requirements

### HG00731 RNA-seq
- **WASP2-processed BAM**: `benchmarking/star_wasp_comparison/results/wasp2rust_fair_snv_*/wasp_filtered.bam`
- **VCF**: `benchmarking/star_wasp_comparison/data/HG00731_het_only_chr.vcf.gz`
- **Sample**: HG00731

### GM12878 ATAC-seq
- **Original BAM**: `/iblm/netapp/data1/aho/atac/Buenrostro2013/merged/GM12878_ATACseq_50k/GM12878_ATACseq_50k_merged.sorted.bam`
- **Remapped BAM**: `benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-15_03-43-24/remap_keep.bam`
- **VCF**: `/iblm/netapp/data1/aho/variants/NA12878.vcf.gz`
- **Sample**: NA12878

### Reference Genome
- **Path**: `/iblm/netapp/data1/external/GRC38/combined/GDC_hg38/bwa_index/GRCh38.d1.vd1.fa` (includes `.dict` required by GATK)

## Tool Versions

- **WASP2-Rust**: Current version (as of Dec 15, 2025)
- **GATK**: Version in WASP2_dev2 conda environment
- **phASER**: Custom installation in `benchmarking/phaser_tool/phaser/`

## Output Format Details

### timing_results.json
```json
{
  "timestamp": "2025-12-16T10:00:00",
  "threads": 8,
  "datasets": {
    "hg00731_rnaseq": {
      "timing": {
        "wasp2_rust_s": 45.2,
        "gatk_s": 1243.5,
        "phaser_s": 567.8
      }
    },
    "gm12878_atacseq": { ... }
  },
  "speedup": {
    "hg00731": {
      "wasp2_vs_gatk": 27.5,
      "wasp2_vs_phaser": 12.6
    }
  }
}
```

### count_comparison.tsv
```
chrom  pos      wasp2_ref  wasp2_alt  gatk_ref  gatk_alt  phaser_ref  phaser_alt
chr1   12345    15         18         15        18        14          19
chr1   67890    22         20         22        20        21          21
```

### original_vs_remapped.tsv
```
chrom  pos    ref_allele  alt_allele  orig_ref  orig_alt  remap_ref  remap_alt  orig_ref_ratio  remap_ref_ratio
chr1   12345  A           G           20        10        15         15         0.667           0.500
chr1   67890  C           T           25        15        20         20         0.625           0.500
```

## Troubleshooting

### Missing Input Files
If the benchmark script reports missing files, check:
1. VCF files are indexed (`.tbi` files present)
2. BAM files are indexed (`.bai` files present)
3. Reference genome is accessible

### phASER Fails
phASER can be finicky. Common issues:
- Requires Python 2.7 (use appropriate conda environment)
- May fail on chromosomes without phased variants
- Output format varies by version

### Slow Bias Comparison
The bias comparison script uses pysam pileup, which can be slow for large BAMs. Consider:
- Running on a subset of chromosomes for testing
- Using only variants with reasonable coverage (>10x)

## Notes

- All paths in scripts use absolute paths from `REPO_ROOT` for SGE compatibility
- Scripts include error handling and will continue if individual tools fail
- Timing uses `date +%s.%N` for sub-second precision
- All plots follow Nature Methods style guidelines (fonts, colors, DPI)

## Citation

If you use these benchmarks, please cite:
- WASP2: [Citation pending]
- GATK: McKenna et al. (2010)
- phASER: Castel et al. (2016)
