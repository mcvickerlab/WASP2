# WASP2 (Currently in development): Allele-specific pipeline for unbiased read mapping, allelic-imbalance analysis, and QTL discovery(WIP)

&nbsp;
## Requirements
- Python >= 3.7
- numpy
- pandas
- polars
- scipy
- pysam
- pybedtools
- typer


## Installation
Recommended installation through conda, and given environment
```shell script
conda env create -f environment.yml
```

&nbsp;
## Allelic Imbalance Analysis
Analysis pipeline currently consists of two tools (Count and Analysis)

&nbsp;
### Count Tool
Process allele specific read counts per SNP.\
Sample names can be provided in order to filter out non-heterozygous SNPs.
ATAC-seq peaks can also be provided to include SNPs that overlap regions of interest.\
Providing samples and regions is highly recommended for allelic-imbalance analysis

**Usage**
```shell script
python WASP2/src/counting count-variants [BAM] [VCF] {OPTIONS}
```

**Required Arguments**
- BAM file containing aligned reads.
- VCF file containing SNP info


**Optional Arguments**
- -s/--samples: Filter SNPs whose genotypes are heterozygous in one or more samples. Accepts comma delimited string, or file with one sample per line. 
- -r/--region: Filter SNPs that overlap peaks/regions of interest. Accepts files in narrowPeak or BED format.
- -o/--out_file: Output file for counts. Defaults to counts.tsv
- -t/--temp_loc: Write intermediary files to a directory instead of deleting. Useful for debugging issues.
- --use_region_names: If regions are provided use region names as identifiers instead of coordinates. Names are denoted in fourth column of BED. Ignored if no name column in BED file.


&nbsp;
### Analysis Tool
Analyzes Allelic Imbalance per ATAC peak given allelic count data

**Usage**
```shell script
python run_analysis.py analysis [COUNTS] {OPTIONS}
```
**Required Arguments**
- COUNTS: first positional argument, output data from count tool

**Optional Arguments**
- --min: Minimum allele count needed for analysis. (Default. 10)
- -o/--output: Directory to output counts. Defaults to CWD if not given. (Default. CWD)
- -m/--model: Model used for measuring imbalance. Choice of "single", "linear", or "binomial". (Default. "single")


&nbsp;
## Unbiased Allele-Specific Read Mapping
Mappability filtering pipeline for correcting allelic mapping biases.\
First, reads are mapped normally using a mapper chosen by the user (output as BAM). Then mapped reads that overlap single nucleotide polymorphisms (SNPs) are identified. For each read that overlaps a SNP, its genotype is swapped with that of the other allele and the read is re-mapped. Re-mapped reads that fail to map to exactly the same location in the genome are discarded.


### Step 1: Create Reads for Remapping
This step identifies reads that overlap snps and creates reads with swapped alleles.

**Usage**
```shell script

python WASP2/src/mapping make-reads [BAM] [VCF] {OPTIONS}
```


**Required Arguments**
- BAM file containing aligned reads.
- VCF file containing SNP info


**Optional Arguments**
- -s/--samples: Filter Polymorphic SNPs in one or more samples. Accepts comma delimited string, or file with one sample per line. 
- -o/--out_dir: Output directory for data to be remapped
- -t/--temp_loc: Write intermediary files to a directory instead of deleting. Useful for debugging issues.
- -j/--out_json: Output json containing wasp file info to this file instead of default. Defaults to [BAM_PREFIX]_wasp_data_files.json


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
python WASP2/src/mapping filter-remapped "${prefix}_remapped.bam" --json "${prefix}_wasp_data_files.json"
```

OR

```shell script
python WASP2/src/mapping filter-remapped "${prefix}_remapped.bam" "${prefix}_to_remap.bam" "${prefix}_keep.bam"
```

**Required Arguments**
- Remapped BAM File
- Either: json or to_remap_bam + keep.bam
    - -j/--json: json containing wasp file info. Default output from make-reads: [BAM_PREFIX]_wasp_data_files.json
    - to_remap_bam: to_remap_bam used to generate swapped alleles. Default: [BAM_PREFIX]_to_remap.bam
    - keep_bam: BAM containing reads that were not remapped. Default: [BAM_PREFIX]_keep.bam

**Optional Arguments**
- -o/--out_bam: File to write filtered bam. Defaults to [BAM_PREFIX]_wasp_filt.bam.
-  --remap_keep_bam: Output bam file with kept reads to this file if provided.
-  --remap_keep_file: Output txt file with kept reads names to this file if provided.


## Future Updates

- Count and Analysis
    - Need to implement RNA-Seq and Gene support 
    - Update Analysis CLI to better work with new counts (Previous analysis CLI in feat-singlecell branch)
    - Reimplement single-cell ssupport and add to CLI parser (Previous Single-Cell counting in feat-singlecell branch)

- Remapping
    - Need to implement single-end and unphased data support
    - Add rmdup step into pipeline
    - More optimization needed

