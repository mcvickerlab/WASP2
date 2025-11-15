# Counting Module - Technical Documentation

**Module**: `src/counting/`
**Purpose**: Extract allele-specific read counts from sequencing data
**Lines of Code**: 1,430 (7 files)
**Generated**: 2025-11-15 (Phase 1.2)

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [File-by-File Analysis](#file-by-file-analysis)
4. [Key Algorithms](#key-algorithms)
5. [Data Structures](#data-structures)
6. [API Reference](#api-reference)
7. [Dependencies](#dependencies)
8. [Known Issues](#known-issues)

---

## Overview

### Purpose

The Counting module processes aligned sequencing reads (BAM format) and variant calls (VCF format) to:
1. Identify reads overlapping heterozygous SNPs
2. Count how many reads support each allele (reference vs alternate)
3. Output allele-specific counts for downstream analysis

### Supported Workflows

1. **Bulk RNA-seq/ATAC-seq**: Count alleles across all reads
2. **Single-cell**: Count alleles per cell barcode (scRNA-seq, scATAC-seq)

### Key Features

- ✅ VCF filtering by sample heterozygosity
- ✅ Region filtering (BED, narrowPeak, GTF, GFF3)
- ✅ Gene annotation integration (exon, transcript, gene levels)
- ✅ Cell barcode tracking (single-cell)
- ✅ Output in TSV or H5AD (AnnData) format

---

## Architecture

### Module Structure

```
src/counting/
│
├── __main__.py (221 lines)
│   └─ CLI Entry Point (Typer)
│       ├─ count_variants()      → Bulk workflow
│       └─ count_variants_sc()   → Single-cell workflow
│
├── BULK WORKFLOW:
│   ├── run_counting.py (230 lines)
│   │   └─ Orchestrator
│   │       ├─ WaspCountFiles (class) - File management
│   │       └─ run_count_variants() - Main pipeline
│   │
│   └── count_alleles.py (123 lines)
│       └─ Core counting logic
│           ├─ make_count_df() - Main entry
│           ├─ count_snp_alleles() - Per-SNP counting
│           └─ find_read_aln_pos() - Binary search helper
│
├── SINGLE-CELL WORKFLOW:
│   ├── run_counting_sc.py (178 lines)
│   │   └─ Orchestrator
│   │       ├─ WaspCountSC (class) - File management
│   │       └─ run_count_variants_sc() - Main pipeline
│   │
│   └── count_alleles_sc.py (228 lines)
│       └─ Core counting logic
│           ├─ make_count_matrix() - AnnData creation
│           ├─ count_bc_snp_alleles() - Per-barcode counting
│           └─ CountStatsSC (class) - Statistics tracking
│
├── SHARED UTILITIES:
│   ├── filter_variant_data.py (237 lines)
│   │   ├─ vcf_to_bed() - VCF → BED conversion
│   │   ├─ intersect_vcf_region() - Region filtering
│   │   ├─ parse_intersect_region_new() - Parse intersections
│   │   └─ gtf_to_bed() - GTF → BED conversion
│   │
│   └── parse_gene_data.py (213 lines)
│       ├─ WaspGeneData (class) - Gene metadata
│       ├─ parse_gene_file() - GTF/GFF3 parsing
│       ├─ make_gene_data() - Create gene BED
│       └─ parse_intersect_genes_new() - Parse gene intersections
│
└── DEPRECATED (Old versions, not called):
    ├─ parse_intersect_region()
    └─ parse_intersect_genes()
```

### Design Pattern: Orchestrator + Worker

**Pattern**:
```
CLI (__main__.py)
  ↓
Orchestrator (run_counting.py)
  ├─ Setup (file paths, validation)
  ├─ VCF filtering (filter_variant_data.py)
  ├─ Gene parsing (parse_gene_data.py) [optional]
  └─ Counting (count_alleles.py)
      ↓
Output (TSV or H5AD)
```

**Benefits**:
- Clear separation of concerns
- Testable components
- Flexible workflow composition

---

## File-by-File Analysis

### 1. `__main__.py` (221 lines)

**Purpose**: CLI interface using Typer framework

**Commands**:

#### 1.1 `count-variants` (Bulk)

```python
@app.command()
def count_variants(
    bam: str,                    # Required: BAM file
    vcf: str,                    # Required: VCF file
    samples: List[str] = None,   # Filter for heterozygous SNPs in samples
    region_file: str = None,     # BED/narrowPeak/GTF/GFF3 regions
    out_file: str = None,        # Output TSV (default: counts.tsv)
    temp_loc: str = None,        # Keep temp files for debugging
    use_region_names: bool = False,  # Use BED name column
    gene_feature: str = None,    # GTF feature type (default: exon)
    gene_attribute: str = None,  # GTF attribute for ID
    gene_parent: str = None      # GTF parent attribute
):
```

**Invocation**:
```bash
python -m src.counting count-variants \
    data.bam variants.vcf \
    --samples NA12878 \
    --region peaks.bed \
    --out counts.tsv
```

#### 1.2 `count-variants-sc` (Single-Cell)

```python
@app.command()
def count_variants_sc(
    bam: str,                    # Required: BAM with CB tags
    vcf: str,                    # Required: VCF file
    barcodes: str,               # Required: Cell barcode file
    samples: List[str] = None,   # Filter SNPs (recommend 1 sample)
    feature_file: str = None,    # BED/narrowPeak features
    out_file: str = None,        # Output H5AD (default: allele_counts.h5ad)
    temp_loc: str = None         # Keep temp files
):
```

**Output**: AnnData H5AD file with:
- `.X`: Total counts (ref + alt + other)
- `.layers['ref']`: Reference allele counts
- `.layers['alt']`: Alternate allele counts
- `.layers['other']`: Other allele counts
- `.obs`: SNP metadata (chrom, pos, ref, alt)
- `.var_names`: Cell barcodes

**Issues Found**:

```python
# Line 118-121: POTENTIAL BUG
if len(samples) > 0:  # TypeError if samples is None!
    samples=samples[0]
else:
    samples=None
```

**Problem**: `len()` on `None` raises TypeError. Should check `if samples is not None and len(samples) > 0`.

**TODOs**:
- Line 16: `# TODO GOTTA TEST THIS` (entire CLI untested?)
- Line 138-139: `# TODO TEST CASES FOR TYPER` / `# TODO UNIT TEST NEXT`
- Line 175: `# TODO: Implement genes gtf/gff format` (for single-cell)

---

### 2. `run_counting.py` (230 lines)

**Purpose**: Orchestrates bulk counting workflow

#### 2.1 `WaspCountFiles` Class

**Purpose**: Manages file paths and parsing logic

```python
class WaspCountFiles:
    def __init__(self, bam_file, vcf_file, region_file=None,
                 samples=None, use_region_names=False,
                 out_file=None, temp_loc=None):
```

**Responsibilities**:
- Parse input file paths
- Determine file types (BED vs GTF vs GFF3)
- Generate intermediate file paths
- Handle sample parsing (file or comma-delimited string)

**File Type Detection** (lines 71-93):
```python
if re.search(r'\.(.*Peak|bed)(?:\.gz)?$', f_ext, re.I):
    self.region_type = "regions"
elif re.search(r'\.g[tf]f(?:\.gz)?$', f_ext, re.I):
    self.region_type = "genes"
    self.is_gene_file = True
elif re.search(r'\.gff3(?:\.gz)?$', f_ext, re.I):
    self.region_type = "genes"
    self.is_gene_file = True
else:
    print("invalid ftype")  # Should raise exception!
```

**Issue**: Uses `print()` instead of raising exception for invalid file types.

**Sample Parsing** (lines 39-48):
```python
if isinstance(self.samples, str):
    if Path(self.samples).is_file():
        # Read from file
        with open(self.samples) as sample_file:
            self.samples = [l.strip() for l in sample_file]
    else:
        # Parse comma-delimited
        self.samples = [s.strip() for s in self.samples.split(",")]
```

**Good**: Flexible input (file or string)
**Issue**: No validation that samples exist in VCF

#### 2.2 `run_count_variants()` Function

**Signature**:
```python
@tempdir_decorator
def run_count_variants(bam_file, vcf_file, region_file=None,
                       samples=None, use_region_names=None,
                       out_file=None, temp_loc=None,
                       gene_feature=None, gene_attribute=None,
                       gene_parent=None):
```

**Pipeline Steps**:

1. **Initialize File Manager** (lines 138-144):
   ```python
   count_files = WaspCountFiles(...)
   ```

2. **Determine Genotype Handling** (lines 147-154):
   ```python
   with_gt = False
   if (count_files.samples is not None) and (len(count_files.samples) == 1):
       with_gt = True  # Include GT column for single sample
   ```

3. **VCF Filtering** (lines 157-161):
   ```python
   vcf_to_bed(vcf_file=count_files.vcf_file,
              out_bed=count_files.vcf_bed,
              samples=count_files.samples,
              include_gt=with_gt)
   ```

4. **Gene Parsing** (lines 169-192, if GTF/GFF3):
   ```python
   if count_files.gtf_bed is not None:
       gene_data = make_gene_data(
           gene_file=count_files.region_file,
           out_bed=count_files.gtf_bed,
           feature=gene_feature,
           attribute=gene_attribute,
           parent_attribute=gene_parent)
       regions_to_intersect = count_files.gtf_bed
   ```

5. **Region Intersection** (lines 190-192):
   ```python
   intersect_vcf_region(vcf_file=count_files.vcf_bed,
                        region_file=regions_to_intersect,
                        out_file=count_files.intersect_file)
   ```

6. **Parse Intersections** (lines 197-219):
   ```python
   if intersect_genes:
       df = parse_intersect_genes_new(...)
   elif with_gt:
       df = parse_intersect_region_new(..., samples=["GT"])
   else:
       df = parse_intersect_region_new(..., samples=None)
   ```

7. **Count Alleles** (lines 224-225):
   ```python
   count_df = make_count_df(bam_file=count_files.bam_file, df=df)
   ```

8. **Write Output** (line 228):
   ```python
   count_df.write_csv(count_files.out_file, has_header=True, separator="\t")
   ```

**Decorator: `tempdir_decorator`** (lines 106-121):

```python
@functools.wraps(func)
def tempdir_wrapper(*args, **kwargs):
    if kwargs.get("temp_loc", None) is not None:
        return func(*args, **kwargs)  # Use provided temp_loc
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            kwargs["temp_loc"] = tmpdir
            return func(*args, **kwargs)  # Auto cleanup
```

**Purpose**: Manages temporary directory lifecycle
**Good**: Automatic cleanup if no temp_loc provided
**Issue**: Doesn't guarantee cleanup on error in all paths

**TODOs**:
- Line 15: `# Should I put this in separate file?` (WaspCountFiles class)
- Line 98: `# TODO UPDATE THIS WHEN I ADD AUTOPARSERS`
- Line 164: `# TODO PARSE GENE FEATURES AND ATTRIBUTES`
- Line 174: `# TODO UPDATE THIS WHEN I ADD AUTOPARSERS AND VALIDATORS`
- Line 196: `# TODO validate`
- Line 221: `# Should I include a filt bam step???`
- Line 230: `# Should i return for use in analysis pipeline?`

---

### 3. `count_alleles.py` (123 lines)

**Purpose**: Core allele counting logic (bulk)

#### 3.1 `find_read_aln_pos()` - Binary Search Helper

```python
def find_read_aln_pos(read, pos):
    aln_list = read.get_aligned_pairs(True)
    i = bisect_left(aln_list, pos, key=lambda x: x[1])

    if i != len(aln_list) and aln_list[i][1] == pos:
        return aln_list[i][0]  # Query position
    else:
        return None
```

**Purpose**: Fast O(log n) lookup of query position for reference position
**Input**: `read` (pysam.AlignedSegment), `pos` (genomic position)
**Output**: Query position (index into read.query_sequence) or None
**Algorithm**: Binary search on aligned pairs

#### 3.2 `make_count_df()` - Main Entry Point

```python
def make_count_df(bam_file, df):
    """
    Make DF containing all intersections and allele counts

    :param str bam_file: Path to BAM file
    :param DataFrame df: Dataframe of intersections from parse_intersect_*()
    :return DataFrame: DataFrame of counts
    """
```

**Process**:

1. **Get Unique Chromosomes** (lines 33-34):
   ```python
   chrom_list = df.get_column("chrom").unique(maintain_order=True)
   ```

2. **Open BAM** (line 38):
   ```python
   with AlignmentFile(bam_file, "rb") as bam:
   ```

3. **Process Each Chromosome** (lines 40-54):
   ```python
   for chrom in chrom_list:
       chrom_df = df.filter(pl.col("chrom") == chrom)
       snp_list = chrom_df.select(["pos", "ref", "alt"]).unique(...)

       try:
           count_list.extend(count_snp_alleles(bam, chrom, snp_list))
       except ValueError:
           print(f"Skipping {chrom}: Contig not found\n")
   ```

4. **Create Count DataFrame** (lines 61-71):
   ```python
   chrom_enum = pl.Enum(df.get_column("chrom").cat.get_categories())
   count_df = pl.DataFrame(
       count_list,
       schema={"chrom": chrom_enum,
               "pos": pl.UInt32,
               "ref_count": pl.UInt16,
               "alt_count": pl.UInt16,
               "other_count": pl.UInt16})
   ```

5. **Join with Input DF** (lines 74-75):
   ```python
   df = df.with_columns([pl.col("chrom").cast(chrom_enum)]
                       ).join(count_df, on=["chrom", "pos"], how="left")
   ```

**Performance**: Uses `timeit` for timing (lines 36, 47, 54, 57-58)

**Good Design**:
- ✅ Processes by chromosome (memory efficient)
- ✅ Uses Polars for fast DataFrame operations
- ✅ Catches `ValueError` for missing contigs

**Issues**:
- ⚠️ Uses `print()` instead of logging
- ⚠️ No progress indicator for large files

#### 3.3 `count_snp_alleles()` - Per-SNP Counting

```python
def count_snp_alleles(bam, chrom, snp_list):
    """Helper function called by make_count_df()"""

    read_set = set()  # Prevent double-counting
    allele_counts = []

    for pos, ref, alt in snp_list:
        ref_count, alt_count, other_count = 0, 0, 0

        # Fetch reads overlapping SNP
        for read in bam.fetch(chrom, pos-1, pos):

            # Skip if already counted
            if read.query_name in read_set:
                continue

            read_set.add(read.query_name)
            seq = read.query_sequence

            # CRITICAL: Linear search (O(n) instead of O(log n)!)
            for qpos, refpos in read.get_aligned_pairs(True):
                if refpos == pos-1:
                    if seq[qpos] == ref:
                        ref_count+=1
                    elif seq[qpos] == alt:
                        alt_count+=1
                    else:
                        other_count+=1
                    break  # Found position, stop searching

        allele_counts.append((chrom, pos, ref_count, alt_count, other_count))

    return allele_counts
```

**🚨 CRITICAL PERFORMANCE BUG** (lines 107-120):

```python
# TODO Update with binary search
for qpos, refpos in read.get_aligned_pairs(True):  # O(n) linear search!
    if refpos == pos-1:
        ...
```

**Problem**: Uses O(n) linear search even though `find_read_aln_pos()` binary search function EXISTS in same file!

**Should be** (like single-cell version does):
```python
qpos = find_read_aln_pos(read, pos-1)  # O(log n)
if qpos is not None:
    if seq[qpos] == ref:
        ref_count += 1
    ...
```

**Impact**: Slower on reads with many indels (long aligned pairs list)

**Evidence**: Single-cell version (count_alleles_sc.py:209) DOES use binary search!

**Other Issues**:

1. **Read Set Behavior** (line 88, 100-103):
   ```python
   read_set = set()  # Outside SNP loop!

   for pos, ref, alt in snp_list:
       ...
       if read.query_name in read_set:
           continue
       read_set.add(read.query_name)
   ```

   **Concern**: `read_set` is shared across ALL SNPs on chromosome. A read with 2 SNPs will only be counted for the FIRST SNP encountered.

   **Is this intentional?** Depends on whether reads should be counted at multiple SNPs or not. Needs clarification.

2. **Commented Code** (line 93):
   ```python
   # read_set = set()  # Commented out - was this meant to be inside loop?
   ```

---

### 4. `filter_variant_data.py` (237 lines)

**Purpose**: VCF filtering and intersection with regions

#### 4.1 `vcf_to_bed()` - VCF Filtering via bcftools

```python
def vcf_to_bed(vcf_file, out_bed, samples=None, include_gt=True):
```

**Process**:

1. **Base VCF View Command**:
   ```bash
   bcftools view {vcf_file} -m2 -M2 -v snps -Ou
   ```
   - `-m2 -M2`: Only biallelic variants
   - `-v snps`: Only SNPs
   - `-Ou`: Uncompressed BCF output

2. **Sample-Specific Filtering**:

   **No samples** (lines 28-34):
   ```bash
   bcftools view --drop-genotypes
   bcftools query -f "%CHROM\t%POS0\t%END\t%REF\t%ALT\n"
   ```

   **Multiple samples** (lines 42-48):
   ```bash
   bcftools view -s {samples} --min-ac 1 --max-ac {2*N-1}
   ```
   - Keeps variants polymorphic in at least one sample

   **Single sample** (lines 50-66):
   ```bash
   bcftools view -s {sample}
   bcftools view --genotype het
   ```
   - Filters for heterozygous SNPs only

3. **Query Output**:
   ```bash
   bcftools query -o {out_bed} -f "{format}"
   ```

   **Without GT**: `%CHROM\t%POS0\t%END\t%REF\t%ALT\n`
   **With GT**: `%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%GT]\n`

**Output Format**: BED-like (0-based start, 1-based end)

**Dependencies**: `bcftools` (subprocess)

**Issue**: No error handling if bcftools fails

#### 4.2 `intersect_vcf_region()` - bedtools Intersection

```python
def intersect_vcf_region(vcf_file, region_file, out_file):
    intersect_cmd = ["bedtools", "intersect",
                     "-a", str(vcf_file),
                     "-b", str(region_file),
                     "-wb"]  # Write both A and B

    with open(out_file, "w") as file:
        subprocess.run(intersect_cmd, stdout=file, check=True)
```

**Dependencies**: `bedtools` (subprocess)

**Output**: Concatenated columns from VCF BED and region BED

#### 4.3 `parse_intersect_region_new()` - Parse Intersections

```python
def parse_intersect_region_new(intersect_file, samples=None,
                                use_region_names=False, region_col=None):
```

**Purpose**: Parse bedtools intersect output into Polars DataFrame

**Logic**:

1. **Base VCF Columns** (lines 130-132):
   ```python
   vcf_cols = ["chrom", "pos0", "pos", "ref", "alt"]
   vcf_schema = [pl.Categorical, pl.UInt32, pl.UInt32,
                 pl.Categorical, pl.Categorical]
   ```

2. **Add Sample Columns** (lines 135-137):
   ```python
   if samples is not None:
       vcf_cols.extend(samples)
       vcf_schema.extend([pl.Categorical] * len(samples))
   ```

3. **Region Column Handling** (lines 155-172):
   - If `use_region_names=True` and BED has name column (4th column):
     ```python
     df = df.rename({df.columns[vcf_ncols+3]: region_col})
     ```
   - Otherwise, concatenate coordinates:
     ```python
     df = df.with_columns(
         pl.concat_str([chr, start, end], separator="_").alias(region_col)
     )
     ```

**Output**: Polars DataFrame with VCF + region columns

#### 4.4 Dead Code

**`parse_intersect_region()` (lines 180-236)**:
- Old version, never called
- Should be removed

---

### 5. `parse_gene_data.py` (213 lines)

**Purpose**: GTF/GFF3 parsing for gene annotations

#### 5.1 `WaspGeneData` Class

```python
class WaspGeneData:
    def __init__(self, gene_file, feature=None,
                 attribute=None, parent_attribute=None):
        self.gene_file = gene_file
        self.feature = feature
        self.attribute = attribute
        self.parent_attribute = parent_attribute
```

**Purpose**: Store gene annotation metadata

#### 5.2 `parse_gene_file()` - GTF/GFF3 Parser

```python
def parse_gene_file(gene_file, feature=None,
                    attribute=None, parent_attribute=None):
```

**Process**:

1. **Read GTF/GFF3** (lines 46-49):
   ```python
   df = pl.read_csv(gene_file, separator="\t",
                    comment_prefix="#",
                    has_header=False,
                    new_columns=gtf_cols)
   ```

2. **Auto-detect Feature** (lines 52-63):
   ```python
   if feature is None:
       feature_list = df.select(pl.col("feature").unique()).to_series()
       if "exon" in feature_list:
           feature = "exon"
       elif "transcript" in feature_list:
           feature = "transcript"
       elif "gene" in feature_list:
           feature = "gene"
   ```

3. **Auto-detect Attribute** (lines 69-80):
   ```python
   if attribute is None:
       if df.get_column("attribute").str.contains(f"{feature}_id").all():
           attribute = f"{feature}_id"
       elif df.get_column("attribute").str.contains("ID").all():
           attribute = "ID"
       elif df.get_column("attribute").str.contains("Name").all():
           attribute = "Name"
   ```

4. **Auto-detect Parent** (lines 86-97):
   ```python
   if parent_attribute is None:
       if df.get_column("attribute").str.contains("Parent").all():
           parent_attribute = "Parent"
       elif df.get_column("attribute").str.contains("transcript_id").all():
           parent_attribute = "transcript_id"
       elif df.get_column("attribute").str.contains("gene_id").all():
           parent_attribute = "gene_id"
   ```

5. **Extract Attributes** (lines 107-114):
   ```python
   attr_regex = fr'{attribute}[=\s]\"?\'?(.*?)\"?\'?;'
   parent_regex = fr'{parent_attribute}[=\s]\"?\'?(.*?)\"?\'?;'

   df = df.with_columns(
       pl.col("start").sub(1),  # Convert to 0-based
       pl.col("attribute").str.extract(attr_regex).alias(attribute),
       pl.col("attribute").str.extract(parent_regex).alias(parent_col)
   ).select(["seqname", "start", "end", attribute, parent_col])
   ```

**Smart Design**: Auto-detects GTF vs GFF3 conventions
**Issue**: No validation if extraction fails (regex returns None)

#### 5.3 `make_gene_data()` - Create BED from GTF

```python
def make_gene_data(gene_file, out_bed, feature=None,
                   attribute=None, parent_attribute=None):
```

**Purpose**: Wrapper that parses GTF and writes BED

**Output BED Format**:
```
seqname  start  end  {attribute}  {parent_attribute}
chr1     1000   2000 exon_001     transcript_001
```

#### 5.4 `parse_intersect_genes_new()` - Parse Gene Intersections

```python
def parse_intersect_genes_new(intersect_file,
                               attribute=None, parent_attribute=None):
```

**Purpose**: Parse bedtools output after intersecting VCF with gene BED

**Expected Columns**:
- 10 columns: VCF(5) + Gene BED(5)
- 11 columns: VCF(6 with GT) + Gene BED(5)

**Output Schema** (without GT):
```python
["chrom", "pos", "ref", "alt", attribute, parent_attribute]
```

**Output Schema** (with GT):
```python
["chrom", "pos", "ref", "alt", "GT", attribute, parent_attribute]
```

**Dead Code**: `parse_intersect_genes()` (old version, lines 156-179)

---

### 6. `run_counting_sc.py` (178 lines)

**Purpose**: Orchestrates single-cell counting workflow

#### 6.1 `WaspCountSC` Class

**MASSIVE CODE DUPLICATION**: ~80% identical to `WaspCountFiles`!

```python
class WaspCountSC:  # vs WaspCountFiles
    def __init__(self, bam_file, vcf_file, barcode_file,  # ← Only difference
                 feature_file, samples=None, ...):
```

**Differences from WaspCountFiles**:
1. Has `barcode_file` parameter
2. Has `feature_file` instead of `region_file`
3. Default output is `.h5ad` instead of `.tsv`
4. Otherwise IDENTICAL code duplication

**Duplication Examples**:

**Sample parsing** (lines 51-60): IDENTICAL
**File extension regex** (lines 86-106): IDENTICAL
**VCF prefix parsing** (lines 73-74): IDENTICAL

**Refactor Opportunity**: Create base class with shared logic

#### 6.2 `run_count_variants_sc()` - Single-Cell Pipeline

```python
@tempdir_decorator
def run_count_variants_sc(bam_file, vcf_file, barcode_file,
                          feature_file=None, samples=None, ...):
```

**Pipeline Steps**:

1. **Initialize** (lines 131-138):
   ```python
   count_files = WaspCountSC(...)
   print(*vars(count_files).items(), sep="\n")  # Debug code in production!
   ```

2. **VCF Filtering** (lines 144-148):
   ```python
   vcf_to_bed(vcf_file=count_files.vcf_file,
              out_bed=count_files.vcf_bed,
              samples=count_files.samples,
              include_gt=True)  # Always True for sc
   ```

3. **Intersect with Features** (lines 150-152):
   ```python
   intersect_vcf_region(vcf_file=count_files.vcf_bed,
                        region_file=count_files.feature_file,
                        out_file=count_files.intersect_file)
   ```

4. **Parse Intersections** (lines 156-161):
   ```python
   df = parse_intersect_region_new(
       intersect_file=count_files.intersect_file,
       samples=count_files.samples,
       use_region_names=use_region_names,
       region_col=None)
   ```

5. **Load Barcodes** (lines 165-167):
   ```python
   with open(count_files.barcode_file, "r") as file:
       bc_dict = {line.rstrip():i for i, line in enumerate(file)}
   ```

6. **Count Alleles** (lines 170-173):
   ```python
   adata = make_count_matrix(bam_file=count_files.bam_file,
                             df=df, bc_dict=bc_dict,
                             include_samples=count_files.samples)
   ```

7. **Write H5AD** (line 176):
   ```python
   adata.write_h5ad(count_files.out_file)
   ```

**Issues**:
- Line 140: `print(*vars(count_files).items(), sep="\n")` - Debug code in production!
- Line 164: `# TODO: handle case where barcode file contains multiple columns`
- Line 177: `# TODO: include output options, (ie MTX, dense?)`
- Line 31: `# TODO: ALSO ACCEPT .h5` (for barcode file)

---

### 7. `count_alleles_sc.py` (228 lines)

**Purpose**: Core single-cell counting logic

#### 7.1 `CountStatsSC` Class - Statistics Tracking

```python
class CountStatsSC:
    def __init__(self):
        self.ref_count = defaultdict(int)
        self.alt_count = defaultdict(int)
        self.other_count = defaultdict(int)

        # Metadata
        self.num_snps = defaultdict(int)
        self.num_barcodes = defaultdict(int)
        self.reads_counted = defaultdict(int)

        # QC stats
        self.reads_skipped_no_barcode = defaultdict(int)
        self.reads_skipped_barcode_no_index = defaultdict(int)
        self.reads_skipped_prev_counted = defaultdict(int)
```

**Purpose**: Track counts and QC metrics per chromosome

**Method**: `stats_to_df()` converts to Pandas DataFrame

#### 7.2 `make_count_matrix()` - Create AnnData

```python
def make_count_matrix(bam_file, df, bc_dict,
                      include_samples=None, include_features=None):
```

**Process**:

1. **Prepare SNP DataFrame** (lines 68-74):
   ```python
   snp_df_cols = ["chrom", "pos", "ref", "alt"]
   if include_samples is not None:
       snp_df_cols.extend(include_samples)

   snp_df = df.select(snp_df_cols).unique(maintain_order=True).with_row_index()
   ```

2. **Count Per Chromosome** (lines 78-100):
   ```python
   with AlignmentFile(bam_file, "rb") as bam:
       for chrom in chrom_list:
           count_bc_snp_alleles(
               bam=bam, bc_dict=bc_dict, chrom=chrom,
               snp_list=chrom_df.select(...).iter_rows(),
               sc_counts=sc_counts)
   ```

3. **Create Sparse Matrices** (lines 104-124):
   ```python
   sparse_ref = csr_matrix(
       (list(sc_counts.ref_count.values()),
        list(zip(*sc_counts.ref_count.keys()))),
       shape=(snp_df.shape[0], len(bc_dict)),
       dtype=np.uint8)
   ```

   **🚨 CRITICAL ISSUE**: Matrix shape is `(SNPs, Cells)` but should be `(Cells, SNPs)`!

   **Standard Convention**:
   - Scanpy: `adata.X` is (cells, genes)
   - SnapATAC2: `adata.X` is (cells, peaks)
   - **Current**: `adata.X` is (SNPs, cells) ← TRANSPOSED!

4. **Create AnnData** (lines 128-135):
   ```python
   adata = ad.AnnData(
       X=sparse_ref+sparse_alt+sparse_other,
       layers={"ref": sparse_ref,
               "alt": sparse_alt,
               "other": sparse_other})
   ```

5. **Annotate** (lines 139-149):
   ```python
   adata.obs = snp_df.to_pandas()  # SNPs as observations (WRONG!)
   adata.var_names = bc_dict.keys()  # Cells as variables (WRONG!)
   ```

   **Should be**:
   ```python
   adata.obs_names = bc_dict.keys()  # Cells as observations
   adata.var = snp_df.to_pandas()    # SNPs as variables
   ```

6. **Add Metadata** (lines 148-173):
   ```python
   if include_samples is not None:
       adata.uns["samples"] = include_samples

   if "region" in df.columns:
       adata.uns["feature"] = df.join(snp_df, ...).select(...).to_pandas()

   adata.uns["count_stats"] = sc_counts.stats_to_df()
   ```

**Issue**: Line 152: `# TODO: Allow for other features besides 'region'`

#### 7.3 `count_bc_snp_alleles()` - Per-Barcode Counting

```python
def count_bc_snp_alleles(bam, bc_dict, chrom, snp_list, sc_counts):
```

**Process**:

```python
read_set = set()  # Prevent double-counting
bc_set = set()    # Track unique barcodes

for idx, pos, ref, alt in snp_list:
    for read in bam.fetch(chrom, pos-1, pos):

        # Skip if already counted
        if read.query_name in read_set:
            sc_counts.reads_skipped_prev_counted[chrom]+=1
            continue

        # Get barcode
        try:
            read_bc = read.get_tag("CB")
        except KeyError:
            sc_counts.reads_skipped_no_barcode[chrom]+=1
            continue

        # Check barcode in index
        if read_bc not in bc_dict:
            sc_counts.reads_skipped_barcode_no_index[chrom]+=1
            continue

        # Find position using binary search (CORRECT!)
        qpos = find_read_aln_pos(read, pos-1)

        try:
            if seq[qpos] == ref:
                sc_counts.ref_count[(idx, bc_dict[read_bc])] += 1
            elif seq[qpos] == alt:
                sc_counts.alt_count[(idx, bc_dict[read_bc])] += 1
            else:
                sc_counts.other_count[(idx, bc_dict[read_bc])] += 1
        except TypeError:
            continue  # qpos is None
        else:
            read_set.add(read.query_name)
            bc_set.add(read_bc)
```

**Good Design**:
- ✅ Uses binary search (`find_read_aln_pos()`) unlike bulk version!
- ✅ Tracks QC metrics
- ✅ Handles missing barcodes gracefully

**Same Issue as Bulk**: `read_set` shared across all SNPs on chromosome

---

## Key Algorithms

### 1. Binary Search for Aligned Position

**Implementation** (`find_read_aln_pos()`):

```python
def find_read_aln_pos(read, pos):
    aln_list = read.get_aligned_pairs(True)
    i = bisect_left(aln_list, pos, key=lambda x: x[1])

    if i != len(aln_list) and aln_list[i][1] == pos:
        return aln_list[i][0]  # Query position
    else:
        return None  # Position not in read
```

**Complexity**: O(log n) where n = number of aligned bases

**Used**: ✅ Single-cell (`count_alleles_sc.py`)
**NOT Used**: ❌ Bulk (`count_alleles.py`) ← Performance bug!

### 2. VCF Filtering Pipeline

**bcftools Commands**:

```bash
# Step 1: Filter biallelic SNPs
bcftools view -m2 -M2 -v snps

# Step 2a: Filter by sample (optional)
bcftools view -s {samples}

# Step 2b: Filter heterozygous (single sample)
bcftools view --genotype het

# Step 3: Convert to BED
bcftools query -f "%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%GT]\n"
```

### 3. GTF Attribute Extraction

**Regex Pattern**:
```python
attr_regex = fr'{attribute}[=\s]\"?\'?(.*?)\"?\'?;'
```

**Matches**:
- GTF: `exon_id "ENSE00001234";`
- GFF3: `ID=exon001;`

**Flexible**: Handles quotes, spaces, semicolons

### 4. Sparse Matrix Construction (Single-Cell)

**Data Structure**:
```python
sc_counts.ref_count = {
    (snp_idx, cell_idx): count,
    ...
}
```

**Conversion to CSR Matrix**:
```python
sparse_ref = csr_matrix(
    (values, (row_indices, col_indices)),
    shape=(n_snps, n_cells),
    dtype=np.uint8
)
```

**Memory Efficient**: Only stores non-zero counts

---

## Data Structures

### Input: Polars DataFrame (Intersections)

**From `parse_intersect_region_new()`**:

```python
# Without regions:
["chrom", "pos", "ref", "alt"]

# With regions:
["chrom", "pos", "ref", "alt", "region"]

# With genes:
["chrom", "pos", "ref", "alt", "exon_id", "transcript_id"]

# With GT:
["chrom", "pos", "ref", "alt", "GT", ...]
```

**Schema**:
```python
{
    "chrom": pl.Categorical,
    "pos": pl.UInt32,
    "ref": pl.Categorical,
    "alt": pl.Categorical,
    "region": pl.Str,  # Optional
    "GT": pl.Categorical,  # Optional
}
```

### Output: Polars DataFrame (Counts)

**Bulk**:
```python
["chrom", "pos", "ref", "alt", "region", "ref_count", "alt_count", "other_count"]
```

**Single-Cell**: AnnData H5AD (see below)

### Output: AnnData (Single-Cell)

```python
adata = ad.AnnData(
    X = sparse_matrix,  # Total counts (ref + alt + other)
    layers = {
        "ref": sparse_ref,
        "alt": sparse_alt,
        "other": sparse_other
    },
    obs = snp_metadata_df,  # ⚠️ Should be cell metadata!
    var_names = cell_barcodes,  # ⚠️ Should be SNP IDs!
    uns = {
        "samples": ["NA12878"],
        "feature": region_snp_mapping,
        "count_stats": qc_stats_df
    }
)
```

**⚠️ CRITICAL ISSUE**: Dimensions are transposed!

---

## Dependencies

### Python Packages

```python
# Data processing
import polars as pl           # DataFrames
import pandas as pd           # Single-cell (AnnData compat)
import numpy as np            # Arrays

# Bioinformatics
from pysam.libcalignmentfile import AlignmentFile  # BAM I/O

# Single-cell
import anndata as ad          # H5AD format
from scipy.sparse import csr_matrix  # Sparse matrices

# System
import tempfile              # Temp directory management
import re                    # Regex for file type detection
from pathlib import Path     # Path operations
from collections import defaultdict, namedtuple
from bisect import bisect_left  # Binary search
```

### System Tools (subprocess)

```bash
bcftools   # VCF filtering (filter_variant_data.py)
bedtools   # Genomic intersections (filter_variant_data.py)
```

**Missing from environment.yml**: Now fixed!

---

## Known Issues

See `docs/modules/COUNTING_ISSUES.md` for complete inventory.

**Critical**:
1. Binary search not used in bulk counting (performance)
2. AnnData dimensions transposed (compatibility)
3. Sample parsing bug (TypeError on None)

**High**:
4. Massive code duplication (WaspCountFiles vs WaspCountSC)
5. Dead code (old parse functions)

**Medium**:
6. Error handling (print instead of raise)
7. Debug code in production
8. 13+ TODO comments

---

**Document Version**: 1.0
**Author**: Phase 1.2 Deep Dive
**Next**: COUNTING_ISSUES.md (Technical Debt Inventory)
