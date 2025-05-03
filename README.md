<h1 align="center">
<img src="doc/wasp2_hex_logo_v1.png" width="300">
</h1>
&nbsp;

# WASP2: Allele-specific pipeline for unbiased read mapping and allelic-imbalance analysis

WASP2 is a comprehensive toolkit designed to perform allele-specific analyses by correcting mapping biases and evaluating allelic imbalance. It supports both bulk and single-cell workflows and is organized into three main command groups:
- **counting:** Generate allele-specific read counts.
- **analysis:** Analyze allelic imbalance from count data.
- **mapping:** Execute the unbiased mapping pipeline.

## Requirements
- Python >= 3.7
- numpy
- pandas
- polars
- scipy
- pysam
- pybedtools
- typer
- anndata
- samtools
- bcftools


## Install from pypi
Recommended installation through conda, and given environment. You can install WASP2 directly via pip:

```shell script
pip install WASP2
```

# Command-Line Interface (CLI)
WASP2’s CLI is organized into three main groups:

## Global Usage
```shell script
WASP2 [OPTIONS] COMMAND [ARGS]...
```

For example, to view the top-level help:
```shell script
WASP2 --help
```
This displays global options (such as --install-completion) and the available command groups: counting, analysis, and mapping.

## Counting Commands
To view help for the counting group:
```shell script
WASP2 counting --help
```

## Analysis Commands
To view help for the analysis group:
```shell script
WASP2 analysis --help
```


&nbsp;
## Allelic Imbalance Analysis
Analysis pipeline currently consists of two tools (Count and Analysis)

&nbsp;
### Count Tool
Process allele specific read counts per SNP.\
Sample names can be provided in order to filter out non-heterozygous SNPs.
Genes and ATAC-seq peaks can also be provided to include SNPs that overlap regions of interest.\
Providing samples and regions is highly recommended for allelic-imbalance analysis

**Usage**
```shell script
python WASP2 counting count-variants [BAM] [VCF] {OPTIONS}
```

**Required Arguments**
- BAM file containing aligned reads.
- VCF file containing SNP info


**Optional Arguments**
- `-s/--samples`: Filter SNPs whose genotypes are heterozygous in one or more samples. Accepts a comma-delimited string, or a file with one sample per line.
- `-r/--region`: Filter SNPs that overlap peaks/regions of interest. Accepts files in narrowPeak, BED, GTF, or GFF3 format.
- `-o/--out_file`: Output file for counts. Defaults to `counts.tsv`.
- `-t/--temp_loc`: Write intermediary files to a directory instead of deleting them. Useful for debugging issues.
- `--use_region_names`: If regions are provided, use region names as identifiers instead of coordinates. Names are denoted in the fourth column of the BED. Ignored if no name column is present.

**RNA-Seq Specific Arguments**
- `--gene_feature`: Feature type in GTF/GFF3 for counting intersecting SNPs. Defaults to `'exon'` for SNP counting.
- `--gene_attribute`: Attribute name from the GTF/GFF3 attribute column to use as ID. Defaults to `'<feature>_id'` in GTF and `'ID'` in GFF3.
- `--gene_parent`: Parent attribute in GTF/GFF3 for the feature used in counting. Defaults to `'transcript_id'` in GTF and `'Parent'` in GFF3.

&nbsp;
### Analysis Tool
Analyzes Allelic Imbalance per ATAC peak given allelic count data

**Usage**
```shell script
python WASP2 analysis find-imbalance [COUNTS] {OPTIONS}
```
**Required Arguments**
- COUNTS: Output data from count tool

**Optional Arguments**
- `-o/--out_file`: Output file to write analysis results to. (Default: `ai_results.tsv`)
- `--min`: Minimum allele count needed for analysis. (Default: `10`)
- `-p/--pseudocount`: Pseudocount added when measuring allelic imbalance. (Default: `1`)
- `--phased`: Calculate allelic imbalance using the phased haplotype model. By default, calculates AI assuming unphased/equal likelihood for each haplotype.
- `--region_col`: Name of the region column for current data. Use `'region'` for ATAC-seq. Plans for `'genes'` for RNA-seq and `'SNP'` for per SNP. Recommended to leave blank. (Default: Auto-parses if none provided)
- `--groupby`: Report allelic imbalance by parent group instead of feature level in RNA-seq counts. Specify the name of the parent column. Not valid if no parent column is present or when using ATAC-seq peaks.  (Default: Report by feature level instead of parent level)

&nbsp;
## Unbiased Allele-Specific Read Mapping
Mappability filtering pipeline for correcting allelic mapping biases.\
First, reads are mapped normally using a mapper chosen by the user (output as BAM). Then mapped reads that overlap single nucleotide polymorphisms (SNPs) are identified. For each read that overlaps a SNP, its genotype is swapped with that of the other allele and the read is re-mapped. Re-mapped reads that fail to map to exactly the same location in the genome are discarded.


### Step 1: Create Reads for Remapping
This step identifies reads that overlap snps and creates reads with swapped alleles.

**Usage**
```shell script

python WASP2 mapping make-reads [BAM] [VCF] {OPTIONS}
```


**Required Arguments**
- BAM file containing aligned reads.
- VCF file containing SNP info


**Optional Arguments**
- `--threads`: Threads to allocate.
- `-s/--samples`: Filter polymorphic SNPs in one or more samples. Accepts a comma-delimited string, or a file with one sample per line.
- `-o/--out_dir`: Output directory for data to be remapped.
- `-t/--temp_loc`: Write intermediary files to a directory instead of deleting them. Useful for debugging issues.
- `-j/--out_json`: Output JSON containing WASP file info to this file instead of the default. Defaults to `[BAM_PREFIX]_wasp_data_files.json`.


### Step 2: Remap Reads
Remap fastq reads using mapping software of choice


**Example**
```shell script
bwa mem -M "BWAIndex/genome.fa" "${prefix}_swapped_alleles_r1.fq" "${prefix}_swapped_alleles_r2.fq" | samtools view -S -b -h -F 4 - > "${prefix}_remapped.bam"
samtools sort -o "${prefix}_remapped.bam" "${prefix}_remapped.bam"
samtools index "${prefix}_remapped.bam"
```


### Step 3: Filter Reads that Fail to Remap
Identify and remove reads that failed to remap to the same position. Creates allelic-unbiased bam file

**Usage**
```shell script
python WASP2 mapping filter-remapped "${prefix}_remapped.bam" --json "${prefix}_wasp_data_files.json"
```

OR

```shell script
python WASP2 mapping filter-remapped "${prefix}_remapped.bam" "${prefix}_to_remap.bam" "${prefix}_keep.bam"
```

**Required Arguments**
- Remapped BAM File
- Either: json or to_remap_bam + keep.bam
- `-j/--json`: JSON containing WASP file info. Default output from make-reads: `[BAM_PREFIX]_wasp_data_files.json`.
- `to_remap_bam`: BAM used to generate swapped alleles. Default: `[BAM_PREFIX]_to_remap.bam`.
- `keep_bam`: BAM containing reads that were not remapped. Default: `[BAM_PREFIX]_keep.bam`.

**Optional Arguments**
- `--threads`: Threads to allocate.
- `-o/--out_bam`: File to write filtered BAM. Defaults to `[BAM_PREFIX]_wasp_filt.bam`.
- `--remap_keep_bam`: Output BAM file with kept reads, if provided.
- `--remap_keep_file`: Output TXT file with kept read names, if provided.

&nbsp;
## Single-Cell Allelic Counts

Process allele specific read counts for single-cell datasets.\
Output counts as anndata containing cell x SNP count matrix.

**Usage**
```shell script
python WASP2 counting count-variants-sc [BAM] [VCF] [BARCODES] {OPTIONS}
```

**Required Arguments**
- `BAM`: BAM file containing aligned reads.
- `VCF`: VCF file containing SNP info.
- `BARCODES`: Barcode file used as index; contains one cell barcode per line.

**Optional Arguments**
- `-s/--samples`: Filter SNPs whose genotypes are heterozygous in one or more samples. Accepts a comma-delimited string, or a file with one sample per line. **RECOMMENDED TO USE ONE SAMPLE AT A TIME.**
- `-f/--feature`: Features used in single-cell experiments. Filter SNPs that overlap regions/features of interest. Accepts BED-formatted files.
- `-o/--out_file`: Output file for counts. Defaults to `allele_counts.h5ad`.
- `-t/--temp_loc`: Write intermediary files to a directory instead of deleting them. Useful for debugging issues.

&nbsp;
## Single-Cell Allelic Imbalance

Estimate allele-specific chromatin acccessibility using single-cell allelic counts.\
Allelic-Imbalance is estimated on a per-celltype basis.

**Usage**
```shell script
python WASP2 counting find-imbalance-sc [COUNTS] [BARCODE_MAP] {OPTIONS}
```

**Required Arguments**
- `COUNTS`: File (.h5ad) containing a matrix of single-cell allelic counts.
- `BARCODE_MAP`: Two-column TSV file mapping specific cell barcodes to a group/celltype. Each line follows the format: `[BARCODE] \t [CELLTYPE]`.

**Optional Arguments**
- `-o/--out_file`: Output file to write analysis results to. (Default: `ai_results_[GROUP].tsv`)
- `--min`: Minimum allele count needed for analysis. (Default: `10`)
- `-p/--pseudocount`: Pseudocount added when measuring allelic imbalance. (Default: `1`)
- `-s/--sample`: Use heterozygous genotypes for this sample in the count matrix. Automatically parsed if data contains 0 or 1 sample.  

**REQUIRED IF MULTIPLE SAMPLES IN DATA.**
- `--phased`: Calculate allelic imbalance using the phased haplotype model. By default, calculates AI assuming unphased/equal likelihood for each haplotype.
- `--unphased`: Explicitly use the unphased model.
- `-z/--z_cutoff`: Remove SNPs and associated regions whose counts exceed the Z-score cutoff.  Provides an extra layer of QC for single-cell allelic counts.

&nbsp;
## Single-Cell Comparative Imbalance

Compare differential allelic-imbalance between celltypes/groups.

**Usage**
```shell script
python WASP2 counting compare-imbalance [COUNTS] [BARCODE_MAP] {OPTIONS}
```

**Required Arguments**
- `COUNTS`: File (.h5ad) containing a matrix of single-cell allelic counts.
- `BARCODE_MAP`: Two-column TSV file mapping specific cell barcodes to a group/celltype. Each line follows the format: `[BARCODE] \t [CELLTYPE]`.

**Optional Arguments**
- `-o/--out_file`: Output file to write analysis results to. (Default: `ai_results_[GROUP1]_[GROUP2].tsv`)
- `--groups/--celltypes`: Specific groups in the barcode map to compare differential allelic imbalance. If provided, requires at least 2 groups; otherwise, compares all group combinations.
- `--min`: Minimum allele count needed for analysis. (Default: `10`)
- `-p/--pseudocount`: Pseudocount added when measuring allelic imbalance. (Default: `1`)
- `-s/--sample`: Use heterozygous genotypes for this sample in the count matrix. Automatically parsed if the data contains 0 or 1 sample.  REQUIRED IF MULTIPLE SAMPLES IN DATA.
- `--phased`: Calculate allelic imbalance using the phased haplotype model. By default, calculates AI assuming unphased/equal likelihood for each haplotype.
- `--unphased`: Explicitly use the unphased model.
- `-z/--z_cutoff`: Remove SNPs and associated regions whose counts exceed the Z-score cutoff. Provides an extra layer of QC for single-cell allelic counts.
