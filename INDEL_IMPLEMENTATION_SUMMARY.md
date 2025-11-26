# WASP2 Indel Support Implementation Summary

## Overview

This document summarizes the implementation of **Stage 1: Variable-Length Indel Support** for WASP2, as outlined in the comprehensive analysis. This implementation enables WASP2 to handle insertions and deletions (indels) in addition to SNPs, using a variable-length approach that is proven, low-risk, and provides immediate benefit.

## Implementation Status

### ✅ Completed Components

1. **CLI Interface** (`src/mapping/__main__.py`)
   - Added `--indels/--snps-only` flag (default: `--snps-only` for backward compatibility)
   - Added `--max-indel-len` parameter (default: 10bp)
   - Added `--insert-qual` parameter (default: Q30)
   - Added `--max-seqs` parameter (default: 64)
   - Added `--same-locus-slop` parameter to `filter-remapped` command (default: 0bp)

2. **Variant Filtering** (`src/wasp2/io/`)
   - Updated `VCFSource.to_bed()` to support indels
   - Updated `compat.py::variants_to_bed()` to pass indel parameters
   - Updated legacy `_vcf_to_bed_bcftools()` function
   - Added bcftools indel length filtering: `strlen(REF)-strlen(ALT)<=max_indel_len`

3. **Parameter Threading** (`src/mapping/`)
   - `intersect_variant_data.py::vcf_to_bed()` - accepts and passes indel params
   - `run_mapping.py::run_make_remap_reads()` - accepts and threads indel params
   - `run_mapping.py::run_wasp_filt()` - accepts and threads same_locus_slop
   - `make_remap_reads.py::write_remap_bam()` - accepts indel params
   - `make_remap_reads.py::swap_chrom_alleles()` - accepts indel params
   - `make_remap_reads.py::swap_chrom_alleles_multi()` - accepts indel params

4. **Position Mapping** (`src/mapping/remap_utils.py`)
   - Implemented `_build_ref2read_maps()` - builds left/right position mappings using `get_aligned_pairs(matches_only=False)`
   - Updated `get_read_het_data()` to support indel-aware position mapping
   - Returns tuple of `(split_seq, split_qual, allele_series)` for indel mode

5. **Quality Score Handling** (`src/mapping/remap_utils.py`)
   - Implemented `_fill_insertion_quals()` - generates quality scores for inserted bases
   - Implemented `make_phased_seqs_with_qual()` - creates haplotypes with quality tracking
   - Updated `write_read()` to accept optional quality scores

### 🔨 Components Requiring Integration

The following components have been updated but need integration work to connect them:

1. **make_remap_reads.py** - Lines 194-237 (swap_chrom_alleles function)
   - Currently uses old SNP-only logic:
     ```python
     r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
     r1_hap_list = [*make_phased_seqs(r1_het_data[0], *r1_het_data[1])]
     ```
   - **Needs update to:**
     ```python
     r1_het_data = get_read_het_data(r1_df, read1, hap_cols,
                                       include_indels=include_indels,
                                       insert_qual=insert_qual)
     if include_indels:
         (r1_hap1_seq, r1_hap1_qual), (r1_hap2_seq, r1_hap2_qual) = make_phased_seqs_with_qual(
             r1_het_data[0], r1_het_data[1], r1_het_data[2][0], r1_het_data[2][1], insert_qual)
         r1_hap_list = [(r1_hap1_seq, r1_hap1_qual), (r1_hap2_seq, r1_hap2_qual)]
     else:
         r1_hap_list = [(seq, None) for seq in make_phased_seqs(r1_het_data[0], *r1_het_data[2])]
     ```

2. **filter_remap_reads.py** - Same-locus slop implementation
   - **Current logic (line ~??):**
     ```python
     if (read1.reference_start, read1.next_reference_start) != pos_dict[read_name]:
         keep_set.remove(read_name)
     ```
   - **Needs update to:**
     ```python
     orig_pos = pos_dict[read_name]
     if abs(read1.reference_start - orig_pos[0]) > same_locus_slop or \
        abs(read1.next_reference_start - orig_pos[1]) > same_locus_slop:
         keep_set.remove(read_name)
     ```

## Key Technical Details

### Position Mapping Algorithm

For indels, we use `get_aligned_pairs(matches_only=False)` which returns `(query_pos, ref_pos)` tuples with `None` for gaps:

- **Matches**: `(5, 100)` - query position 5 aligns to reference position 100
- **Insertions**: `(5, None)` - query position 5 has no reference position (inserted base)
- **Deletions**: `(None, 100)` - reference position 100 has no query position (deleted base)

We build two mappings:
- `ref2q_left`: Maps ref position to nearest left query position
- `ref2q_right`: Maps ref position to nearest right query position

This allows us to handle deletions where the reference position doesn't map to a specific query position.

### Quality Score Handling

For insertions where the alternate allele is longer than the original:
1. Calculate flanking quality scores from left and right segments
2. Average the flanking scores
3. Fill inserted bases with the averaged quality (or `insert_qual` if no flanks available)

For deletions, simply truncate the quality array to match the shorter allele.

### Backward Compatibility

All indel features are **opt-in** via the `--indels` flag:
- Default behavior: `--snps-only` (existing SNP-only functionality preserved)
- Indel mode: `--indels` (enables new indel-aware code paths)

This ensures existing WASP2 users see no breaking changes.

## Usage Examples

### Basic Indel Support

```bash
# Generate remapping reads with indel support
wasp2-map make-reads sample.bam variants.vcf.gz \
  --samples sample1 \
  --indels \
  --max-indel-len 10 \
  --insert-qual 30 \
  --max-seqs 64

# Remap using your favorite aligner
bwa mem genome.fa swapped_alleles_r1.fq swapped_alleles_r2.fq | \
  samtools view -Sb - > remapped.bam
samtools sort -o remapped.sorted.bam remapped.bam
samtools index remapped.sorted.bam

# Filter with same-locus slop for indel micro-homology
wasp2-map filter-remapped remapped.sorted.bam \
  --json sample_wasp_data_files.json \
  --same-locus-slop 2
```

### SNP-Only Mode (Default)

```bash
# Backward compatible - works exactly as before
wasp2-map make-reads sample.bam variants.vcf.gz \
  --samples sample1
```

## Performance Expectations

Based on the analysis:

- **Overhead**: 1.5-2x compared to SNP-only mode (due to quality score handling)
- **Additional reads**: 13-28% more reads retained vs SNP-only (from GTEx data)
- **Memory**: Minimal increase (~5-10% for quality arrays)
- **Disk I/O**: Slightly higher due to longer reads (variable length)

## Parameters Explained

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--indels` | False | Enable indel support (opt-in) |
| `--max-indel-len` | 10 | Maximum indel size (bp) to process |
| `--insert-qual` | 30 | Phred quality score for inserted bases |
| `--max-seqs` | 64 | Max alternate sequences per read (prevents explosion) |
| `--same-locus-slop` | 0 | Tolerance (bp) for same-locus test in filtering |

## Remaining Work

1. **Integration** - Connect the new indel-aware functions in `make_remap_reads.py` (see "Components Requiring Integration" above)
2. **filter_remap_reads.py** - Implement same-locus slop tolerance
3. **Testing** - Validate with synthetic and real data
4. **Documentation** - Update user-facing docs with indel examples

## Validation Plan

### Synthetic Data Tests
1. Create test VCF with known indels (1-10bp insertions/deletions)
2. Create synthetic BAM with reads overlapping indels
3. Run WASP2 with `--indels` and verify:
   - Alternate reads generated correctly
   - Quality scores preserved/filled appropriately
   - Same-locus filtering works with slop

### Real Data Tests
1. Use NA12878 ATAC-seq data (published)
2. Compare indel-aware vs SNP-only modes:
   - Number of reads retained
   - Concordance with known variants
   - Allelic imbalance estimates

### Performance Benchmarks
1. Measure runtime overhead (expect 1.5-2x)
2. Memory profiling
3. Compare to wasp_indel_v1.0 (should be comparable or faster)

## Future Enhancements (Stage 2+)

- **Balanced-trim enumeration** for small indels (1-3bp) in unique regions
- **Repetitive region detection** to skip balanced-trim where it collapses
- **ML-based quality inference** for insertions (better than simple averaging)
- **Multi-threading** for per-chromosome processing

## References

- Original analysis document: "Deep Analysis: Indel Support for WASP2"
- Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12047541/#SD8
- wasp_indel_v1.0: https://github.com/adam-rabinowitz/wasp_indel_v1.0
