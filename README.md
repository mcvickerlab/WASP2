# WASP2 (Currently in pre-development): Allele-specific pipeline for unbiased read mapping(WIP), QTL discovery(WIP), and allelic-imbalance analysis

&nbsp;
## Requirements
- Python >= 3.7
- numpy
- pandas
- scipy
- pysam
- pybedtools

&nbsp;
## Installation
Recommended installation through conda, and given environment
```shell script
conda env create -f environment.yml
```

&nbsp;
## Allelic Imbalance Analysis
Analysis pipeline currently consists of two tools (Count and Analysis)

&nbsp;
## Count Tool
Counts alleles in peaks or genes that overlap heterozygous SNP's

**Usage**
```shell script
python run_analysis.py count -a [BAM] -g [VCF] -s [VCF Sample] -r [Peaks] {OPTIONS}
```

**Required Arguments**
- -a/--alignment: BAM file containing alignments.
- -g/--genotypes: VCF file with genotypes.
- -s/--sample: Sample name in VCF file.
- -r/--regions: Regions of interest
    - ATAC: narrowPeak & BED formatted peaks
    - RNA: GTF formatted gene annotations

**Single-Cell Additional Requirements**
- -sc/--singlecell: Flag that denotes data is single-cell.
- -b/--barcodes: 2 Column TSV that contains barcodes and their group/cell mapping. 

**Optional Arguments**
- --rna/--atac: Denotes if analyzing rna-seq or atac-seq, otherwise infer using inputs
- -ft/--features: (RNA ONLY): NAme of features to analyze in GTF file (eg. transcript, exon, CDS). By default analyzes transcript. Analyzes all features if flag denoted without feature.
- -o/--output: Directory to output counts. (Default. CWD)
- --nofilt: Skip step that pre-filters reads that overlap regions of interest 
- --keeptemps: Keep intermediary files during preprocessing step, outputs to directory if given with flag, otherwise outputs to same location as final output.


&nbsp;
## Analysis Tool
Analyzes Allelic Imbalance per ATAC peak given allelic count data

**Usage**
```shell script
python run_analysis.py analysis [COUNTS] {OPTIONS}
```
**Required Arguments**
- COUNTS: first positional argument, output data from count tool

**Single-Cell Additional Requirements**
- -sc/--singlecell: Flag that denotes data is single-cell

**Optional Arguments**
- --rna/--atac: Denotes if analyzing rna-seq or atac-seq
- -ft/--features: (RNA ONLY): Features to analyze in gtf. By default analyzes all features found separately.
- --min: Minimum allele count needed for analysis. (Default. 10)
- -o/--output: Directory to output counts. Defaults to CWD if not given. (Default. CWD)
- -m/--model: Model used for measuring imbalance. Choice of "single", "linear", or "binomial". (Default. "single")


&nbsp;
## Comparison Tool (EXPERIMENTAL)
Compares Allelic Imbalance between single-cell cells/clusters

**Usage**
```shell script
python run_analysis.py compare [COUNTS] {OPTIONS}
```
**Required Arguments**
- COUNTS: first positional argument, 1 or 2 output files from count tool. 
    - One file will compare cells/clusters in the sample
    - Two files will compare shared cells/clusters between samples

**Optional Arguments**
- --rna/--atac: Denotes if analyzing rna-seq or atac-seq
- --min: Minimum allele count needed for analysis. (Default. 10)
- -o/--output: Directory to output counts. Defaults to CWD if not given. (Default. CWD)


&nbsp;
## TODO
- Unbiased Read Mapping Curently in development


Allelic Imbalance Pipeline
- Counts
    - More robust for different inputs for bulk and single-cell data

- Analysis
    - More specific implementations for single-cell data

- Comparisdon
    - Get Bulk rna/atac and single-cell rna inputs working
    - Allow for analysis with different models
    - Improved stability