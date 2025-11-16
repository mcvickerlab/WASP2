# Mapping Module - Technical Overview

**Phase 1.4 Deep Dive**
**Total LOC:** 1,569 (27.2% of codebase)
**Files:** 7 Python files
**Purpose:** WASP mapping bias removal through allele swapping and remapping

---

## WASP Algorithm Overview

**Problem:** Reference genome mapping bias causes reads with alternative alleles to map less efficiently than reads with reference alleles, confounding allele-specific analysis.

**Solution:**
1. Find reads overlapping heterozygous SNPs
2. Create synthetic reads with alleles swapped to alternative haplotypes
3. Remap synthetic reads with external aligner (e.g., BWA, STAR)
4. Keep only reads that remap to the **same genomic location**
5. Filter removes mapping bias while retaining true signal

---

## Module Structure

### CLI & Orchestration
- **__main__.py** (180 LOC) - 2 commands: `make_reads`, `filter_remapped`
- **run_mapping.py** (240 LOC) - Pipeline orchestrators with decorators
- **wasp_data_files.py** (111 LOC) - File path management class

### Core Processing
- **make_remap_reads.py** (499 LOC) - Allele swapping engine, ReadStats tracking
- **remap_utils.py** (136 LOC) - Sequence manipulation utilities
- **filter_remap_reads.py** (97 LOC) - Post-remapping filter
- **intersect_variant_data.py** (306 LOC) - VCF/BAM intersection creation

---

## Two-Step Pipeline

### Step 1: make_reads
**Input:** BAM + VCF → **Output:** FASTQ files with swapped alleles

```bash
wasp2 mapping make_reads \
    aligned.bam \
    variants.vcf.gz \
    --samples SAMPLE_ID \
    --out_dir ./wasp_output
```

**Internal Flow:**
1. `vcf_to_bed()` - Filter VCF for het SNPs using bcftools
2. `process_bam()` - Extract reads overlapping SNPs
3. `intersect_reads()` - Create intersection BED using pybedtools
4. `swap_chrom_alleles()` - Generate synthetic reads per chromosome
5. Write FASTQ files: `*_swapped_alleles_r1.fq`, `*_swapped_alleles_r2.fq`
6. Export `*_wasp_data_files.json` for step 2

**Allele Swapping:**
```python
# Original read sequence with ref allele "A"
r1_seq = "GCTAGCATGCTAGC"
       #      ^ pos 6 = A (ref)

# Genotype: 0|1 (haplotype 1 has "A", haplotype 2 has "G")
# Create swapped read with alt allele "G"
r1_swapped = "GCTAGCGTGCTAGC"
          #       ^ pos 6 = G (alt)

# New read name encodes original position
new_name = f"{orig_name}_WASP_{r1_pos}_{r2_pos}_{read_num}_{total_reads}"
```

### Step 2: filter_remapped
**Input:** Remapped BAM + wasp_data_files.json → **Output:** WASP-filtered BAM

```bash
# User remaps FASTQ with external aligner (BWA/STAR/etc.)
bwa mem genome.fa swapped_alleles_r1.fq swapped_alleles_r2.fq > remapped.bam

# Filter for reads that remapped to same location
wasp2 mapping filter_remapped \
    remapped.bam \
    --json *_wasp_data_files.json \
    --out wasp_filtered.bam
```

**Internal Flow:**
1. `filt_remapped_reads()` - Check if remapped position == original position
2. `merge_filt_bam()` - Merge filtered remapped reads + keep reads
3. Sort and index final BAM

---

## Key Classes

### WaspDataFiles
Manages all intermediate and output file paths:
```python
class WaspDataFiles:
    bam_file, vcf_file, is_paired, samples, is_phased
    out_dir, temp_loc

    # Intermediate files
    vcf_bed, remap_reads, intersect_file

    # Output files
    to_remap_bam, keep_bam
    remap_fq1, remap_fq2

    def write_data(out_file) # Export to JSON
```

### ReadStats
Tracks reads processed and discarded:
```python
class ReadStats:
    discard_unmapped, discard_improper_pair
    discard_secondary, discard_supplementary
    discard_excess_reads, discard_missing_pair
    remap_pair, write_pair
```

---

## Critical Functions

### swap_chrom_alleles() - Single Sample
Per-chromosome allele swapping for one sample (phased genotypes):
```python
def swap_chrom_alleles(bam_file, out_dir, df, chrom, read_stats):
    # 1. Partition intersection DF by read name
    r1_het_dict = chrom_df.filter(pl.col("mate") == 1).partition_by("read")
    r2_het_dict = chrom_df.filter(pl.col("mate") == 2).partition_by("read")

    # 2. Iterate through read pairs
    for read1, read2 in paired_read_gen_stat(bam, read_stats, chrom):
        # 3. Get heterozygous data per read
        r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
        r2_het_data = get_read_het_data(r2_df, read2, hap_cols)

        # 4. Create phased haplotype sequences
        r1_hap_list = make_phased_seqs(*r1_het_data)  # [hap1_seq, hap2_seq]
        r2_hap_list = make_phased_seqs(*r2_het_data)

        # 5. Write pairs with swapped alleles
        for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list):
            if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq):
                new_name = f"{og_name}_WASP_{r1_pos}_{r2_pos}_{num}_{total}"
                write_read(out_file, read1, r1_hap_seq, new_name)
                write_read(out_file, read2, r2_hap_seq, new_name)

    # 6. Convert BAM → FASTQ using samtools
    subprocess.run(["samtools", "collate", "-u", "-O", out_bam], stdout=PIPE)
    subprocess.run(["samtools", "fastq", "-1", r1_out, "-2", r2_out], input=...)
```

### swap_chrom_alleles_multi() - Multiple Samples
Handles multiple samples with unique haplotype combinations:
```python
def swap_chrom_alleles_multi(bam_file, out_dir, df, chrom, read_stats):
    # Find unique haplotype combinations across all samples
    unique_cols = (
        read_df.select(pl.col(hap_cols).str.concat(""))
        .transpose()
        .unique(subset=["hap"])
        .get_column("column")
    )

    # Generate sequences for each unique haplotype
    r1_hap_list = make_multi_seqs(*r1_het_data)
    r2_hap_list = make_multi_seqs(*r2_het_data)
```

### filt_remapped_reads()
Filters remapped reads by position equality:
```python
def filt_remapped_reads(to_remap_bam, remapped_bam, filt_out_bam):
    pos_dict = {}  # {read_name: (r1_pos, r2_pos)}
    keep_set = set()

    for read1, read2 in paired_read_gen(remapped_bam):
        read_name, r1_pos, r2_pos, num, total = parse_wasp_name(read1.query_name)

        if read_name not in pos_dict:
            pos_dict[read_name] = (r1_pos, r2_pos)
            keep_set.add(read_name)

        # Check position equality
        if (read1.reference_start, read1.next_reference_start) != pos_dict[read_name]:
            keep_set.remove(read_name)  # Failed WASP filter

    # Write kept reads to file
    pysam.view("-N", keep_read_file, "-o", filt_out_bam, to_remap_bam)
```

---

## Decorator Patterns

### @tempdir_decorator
Manages temporary directory creation:
```python
@functools.wraps(func)
def tempdir_wrapper(*args, **kwargs):
    if kwargs.get("temp_loc", None) is not None:
        return func(*args, **kwargs)
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            kwargs["temp_loc"] = tmpdir
            return func(*args, **kwargs)
```

### @check_filt_input
Validates and parses filter_remapped inputs:
```python
@functools.wraps(func)
def filt_wrapper(*args, **kwargs):
    # Accept either JSON or explicit BAM paths
    if kwargs.get("wasp_data_json") is not None:
        json_dict = json.load(wasp_data_json)
        kwargs["to_remap_bam"] = json_dict["to_remap_bam"]
        kwargs["keep_bam"] = json_dict["keep_bam"]
    elif not (kwargs.get("to_remap_bam") and kwargs.get("keep_bam")):
        raise ValueError("Must provide JSON or both BAMs")
```

---

## Data Flow

```
Input BAM + VCF
     ↓
vcf_to_bed (bcftools) → het_snps.bed
     ↓
process_bam (pysam) → to_remap.bam, keep.bam
     ↓
intersect_reads (pybedtools) → intersect.bed
     ↓
make_intersect_df (polars) → DataFrame with genotypes
     ↓
swap_chrom_alleles → swapped_alleles_{chr}.bam (per chrom)
     ↓
samtools collate + fastq → swapped_alleles_r1.fq, r2.fq
     ↓
[USER REMAPS WITH EXTERNAL ALIGNER]
     ↓
filter_remapped_reads → wasp_remap_filt.bam
     ↓
merge_filt_bam → wasp_filtered.bam (final)
```

---

## External Dependencies

| Tool | Usage | Command |
|------|-------|---------|
| bcftools | VCF filtering for het SNPs | `bcftools view/query` |
| pybedtools | BAM/VCF intersections | `BedTool.intersect()` |
| pysam | BAM I/O, view, sort, index, merge | Wrapper around samtools |
| samtools | Collate BAM, convert to FASTQ | `subprocess.run()` |
| polars | Fast DataFrame operations | Lazy evaluation with `.scan_csv()` |

---

## Complexity Metrics

| File | LOC | Functions | Classes | Complexity |
|------|-----|-----------|---------|------------|
| make_remap_reads.py | 499 | 3 | 1 | Very High |
| intersect_variant_data.py | 306 | 5 | 0 | High |
| run_mapping.py | 240 | 4 | 0 | Medium |
| __main__.py | 180 | 2 | 0 | Low |
| remap_utils.py | 136 | 7 | 0 | Medium |
| wasp_data_files.py | 111 | 2 | 1 | Low |
| filter_remap_reads.py | 97 | 2 | 0 | Low |

---

## Current Limitations

**Documented in code:**
- **Single-end reads NOT supported** (line 76: `raise ValueError("Single-End not Implemented")`)
- **Unphased genotypes NOT supported** (line 79: `raise ValueError("Unphased not Implemented")`)
- **Requires sample ID** (line 82: `raise ValueError("Zero samples not supported yet")`)

**Requirements:**
- Paired-end reads only
- Phased VCF (e.g., from SHAPEIT, BEAGLE)
- At least one sample specified

---

## See Also

- **MAPPING_ISSUES.md** - Technical debt and code quality issues
- **COUNTING_MODULE.md** - Downstream allele counting
- **ARCHITECTURE.md** - Overall system design
