# WASP2 Paper Figures

## Directory Structure

```
paper/
├── config.py              # Shared paths and settings
├── run_all.sh             # Master script to regenerate all figures
├── figure1/               # Read Mapping
│   ├── scripts/
│   ├── data/
│   └── plots/
├── figure2/               # Counting
│   ├── scripts/
│   ├── data/
│   └── plots/
├── figure3/               # Statistics
│   ├── scripts/
│   ├── data/
│   └── plots/
├── figure4/               # Single Cell
│   ├── scripts/
│   ├── data/
│   └── plots/
├── supplementary/
└── shared_data/           # Large datasets symlinked here
```

## Figure Descriptions

### Figure 1: Read Mapping
- **1A**: WASP2 pipeline schematic
- **1B**: Speed comparison (WASP2 vs WASP1, ATAC + RNA)
- **1C**: Reads overlapping SNVs vs INDELs (pre/post remapping)

### Figure 2: Counting
- **2A**: Speed comparison (WASP2 vs GATK vs phASER)
- **2B**: Ref/Alt count comparison across methods
- **2C**: Original vs remapped BAM bias analysis

### Figure 3: Statistics
- **3A**: Statistical methods schematic
- **3B**: QQ plots for different methods
- **3C**: SNV-based volcano plot
- **3D**: Gene imprinting analysis (SNVs + INDELs)
- **3E**: Het/homo QTL stratification

### Figure 4: Single Cell
- Single cell allelic imbalance analysis

## Usage

```bash
# Run all figure generation
./run_all.sh

# Run individual figures
python figure1/scripts/generate_figure1.py
python figure2/scripts/generate_figure2.py
python figure3/scripts/generate_figure3.py
```

## Data Sources

| Dataset | Location | Used For |
|---------|----------|----------|
| GM12878 ATAC-seq | `/iblm/netapp/data1/aho/atac/...` | Fig 1, 2, 3 |
| HG00731 RNA-seq | `benchmarking/star_wasp_comparison/` | Fig 1B, 2A |
| NA12878 VCF | `/iblm/netapp/data1/aho/variants/` | All figures |
| iPSCORE QTLs | TBD | Fig 3B, 3E |
