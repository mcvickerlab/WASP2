# WASP2 Full Pipeline Execution Plan

**Created**: 2025-11-15
**Purpose**: Document end-to-end pipeline execution for baseline establishment
**Test Data**: Aaron Ho's WASP test data (NA12878 chr10, ATAC-seq)

---

## 🚨 **CRITICAL FINDINGS - Environment Setup Issues**

### **Missing System Dependencies**

❌ **NOT in environment.yml but REQUIRED by code:**

1. **bcftools** - Required by `src/counting/filter_variant_data.py` (lines 19, 23, 57)
   - Used for VCF filtering via subprocess
   - Called 3+ times in VCF processing pipeline

2. **bedtools** - Required by `src/counting/filter_variant_data.py` (line 112)
   - Used for genomic interval intersections
   - Critical for region filtering

3. **samtools** - Required by `src/mapping/run_mapping.py`
   - Used for BAM sorting and indexing
   - Required for mapping pipeline

### **Environment Status**

```bash
# Current environment check:
bcftools: NOT FOUND ❌
bedtools: NOT FOUND ❌
samtools: NOT FOUND ❌

# Python packages: NOT INSTALLED ❌
# (No conda environment activated)
```

### **Fix Required**

Update `environment.yml`:
```yaml
dependencies:
  - python=3.9.*
  - numpy
  - pandas
  - polars
  - scipy
  - pysam
  - pybedtools
  - bedtools      # ← ADD THIS
  - typer
  - anndata
  # ADD THESE:
  - bcftools      # ← CRITICAL - used heavily
  - samtools      # ← CRITICAL - used by mapping
```

---

## 📊 **Test Data Inventory**

### **Files Available**

| File | Size | Description |
|------|------|-------------|
| `filter_chr10.vcf` | 12 MB | NA12878 chr10 SNPs (111,454 variants) |
| `NA12878_snps_chr10.bed` | 2.6 MB | Matching genomic intervals |
| `CD4_ATACseq_Day1_merged_filtered.sort.bam` | 7.6 MB | Small ATAC-seq BAM subset |
| `CD4_ATACseq_Day1_merged_filtered.sort.bam.bai` | 2.0 MB | BAM index file |
| `as_counts.txt` | 274 bytes | Expected output format example |

### **Data Characteristics**

**VCF File** (`filter_chr10.vcf`):
- Format: VCFv4.1
- Source: Illumina Platinum Genomes (2016-1.0)
- Reference: hg38
- Sample: NA12878
- Variants: 111,454 chr10 SNPs
- File date: 2016-06-03

**BAM File**:
- Format: BGZF (blocked gzip)
- Type: ATAC-seq (CD4 T cells, Day 1)
- Status: Sorted and indexed
- Size: 7.6 MB (compressed)

**Expected Output** (`as_counts.txt`):
```
chrom   pos      ref  alt  peak                    ref_count  alt_count  other_count
chr1    1019397  C    A    chr1_1019383_1019826    0          2          0
chr1    1019636  A    G    chr1_1019383_1019826    0          1          0
```

---

## 🔄 **Full Pipeline Execution Plan**

### **Prerequisites**

1. **Install Conda Environment**:
   ```bash
   cd /home/user/WASP2-exp
   conda env create -f environment.yml
   conda activate WASP2
   ```

2. **Verify Installation**:
   ```bash
   # Check system tools
   bcftools --version
   bedtools --version
   samtools --version

   # Check Python imports
   python -c "import pysam, polars, scipy, typer, anndata; print('All imports OK')"
   ```

3. **Prepare Test Data** (already extracted):
   ```bash
   ls -lh test_data/
   # Should show BAM, VCF, BED files
   ```

---

## 📝 **Pipeline Steps**

### **OPTION 1: Counting Only (Skip Mapping)**

**Assumption**: Test BAM is already WASP-filtered or we're testing counting in isolation

#### **Step 1: Count Alleles**

```bash
python -m src.counting count-variants \
    test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam \
    test_data/filter_chr10.vcf \
    --samples NA12878 \
    --region test_data/NA12878_snps_chr10.bed \
    --out baselines/counting/counts.tsv \
    --temp baselines/counting/temp

# Expected outputs:
# - baselines/counting/counts.tsv (main output)
# - baselines/counting/temp/*.bed (intermediates)
```

**Validation**:
```bash
# Check output format
head baselines/counting/counts.tsv

# Compare to expected format
diff <(head baselines/counting/counts.tsv) <(head test_data/as_counts.txt)

# Record stats
wc -l baselines/counting/counts.tsv
md5sum baselines/counting/counts.tsv > baselines/counting/counts.md5
```

#### **Step 2: Analyze Allelic Imbalance**

```bash
python -m src.analysis find-imbalance \
    baselines/counting/counts.tsv \
    --out baselines/analysis/ai_results.tsv \
    --min 10 \
    --pseudocount 1

# Expected outputs:
# - baselines/analysis/ai_results.tsv (AI statistics)
```

**Validation**:
```bash
# Check output
head baselines/analysis/ai_results.tsv
wc -l baselines/analysis/ai_results.tsv

# Record baseline
md5sum baselines/analysis/ai_results.tsv > baselines/analysis/ai_results.md5

# Check for regions with significant imbalance
awk '$6 < 0.05' baselines/analysis/ai_results.tsv | wc -l
# (assuming column 6 is p-value or FDR)
```

---

### **OPTION 2: Full Pipeline (Including Mapping)**

#### **Step 1: Create Reads for Remapping**

```bash
python -m src.mapping make-reads \
    test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam \
    test_data/filter_chr10.vcf \
    --samples NA12878 \
    --out_dir baselines/mapping \
    --temp baselines/mapping/temp

# Expected outputs:
# - baselines/mapping/*_swapped_alleles_r1.fq
# - baselines/mapping/*_swapped_alleles_r2.fq
# - baselines/mapping/*_to_remap.bam
# - baselines/mapping/*_keep.bam
# - baselines/mapping/*_wasp_data_files.json
```

**Validation**:
```bash
# Check FASTQs were created
ls -lh baselines/mapping/*.fq

# Count reads in each file
echo "=== Swapped reads (FASTQ) ==="
grep -c "^@" baselines/mapping/*_r1.fq

echo "=== Reads to remap (BAM) ==="
samtools view -c baselines/mapping/*_to_remap.bam

echo "=== Reads to keep (BAM) ==="
samtools view -c baselines/mapping/*_keep.bam

# Inspect metadata JSON
cat baselines/mapping/*_wasp_data_files.json
```

#### **Step 2: Remap Swapped Reads**

⚠️ **REQUIRES: Reference genome (hg38) and aligner (BWA recommended)**

```bash
# Example with BWA (requires hg38 reference genome)
REF_GENOME=/path/to/hg38.fa
PREFIX=baselines/mapping/CD4_ATACseq_Day1_merged_filtered

bwa mem -t 4 \
    $REF_GENOME \
    ${PREFIX}_swapped_alleles_r1.fq \
    ${PREFIX}_swapped_alleles_r2.fq \
    | samtools view -b -F 4 - \
    > ${PREFIX}_remapped.bam

# Sort and index
samtools sort -o ${PREFIX}_remapped.sorted.bam ${PREFIX}_remapped.bam
samtools index ${PREFIX}_remapped.sorted.bam
```

**Note**: This step is **blocked** without:
- hg38 reference genome (~3 GB download)
- BWA aligner installed
- Sufficient compute resources

**Alternative**: Skip for now, proceed with existing BAM for counting tests

#### **Step 3: Filter Remapped Reads**

```bash
python -m src.mapping filter-remapped \
    baselines/mapping/CD4_ATACseq_Day1_merged_filtered_remapped.sorted.bam \
    --json baselines/mapping/CD4_ATACseq_Day1_merged_filtered_wasp_data_files.json \
    --out baselines/mapping/wasp_filt.bam

# Expected output:
# - baselines/mapping/wasp_filt.bam (unbiased BAM)
```

**Validation**:
```bash
# Count filtered reads
samtools view -c baselines/mapping/wasp_filt.bam

# Compare to original
echo "Original BAM:"
samtools view -c test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam

echo "After WASP filtering:"
samtools view -c baselines/mapping/wasp_filt.bam

# Calculate filtering rate
```

#### **Step 4: Count with WASP-Filtered BAM**

```bash
python -m src.counting count-variants \
    baselines/mapping/wasp_filt.bam \
    test_data/filter_chr10.vcf \
    --samples NA12878 \
    --region test_data/NA12878_snps_chr10.bed \
    --out baselines/counting/counts_wasp_filtered.tsv
```

#### **Step 5: Analyze**

```bash
python -m src.analysis find-imbalance \
    baselines/counting/counts_wasp_filtered.tsv \
    --out baselines/analysis/ai_results_wasp_filtered.tsv
```

---

## 🎯 **Baseline Artifacts to Save**

### **After Successful Run, Save:**

1. **Mapping Outputs**:
   ```bash
   baselines/mapping/
   ├── *_swapped_alleles_r1.fq (MD5, line count)
   ├── *_swapped_alleles_r2.fq (MD5, line count)
   ├── *_to_remap.bam (MD5, read count)
   ├── *_keep.bam (MD5, read count)
   ├── *_wasp_data_files.json (save entire file)
   └── wasp_filt.bam (MD5, read count)
   ```

2. **Counting Outputs**:
   ```bash
   baselines/counting/
   ├── counts.tsv (MD5, row count, first 10 lines)
   └── counts_wasp_filtered.tsv (if mapping was run)
   ```

3. **Analysis Outputs**:
   ```bash
   baselines/analysis/
   ├── ai_results.tsv (MD5, row count, summary stats)
   └── ai_results_wasp_filtered.tsv (if mapping was run)
   ```

4. **Metadata**:
   ```bash
   baselines/
   ├── execution_log.txt (full command output)
   ├── timing.txt (execution times for each step)
   └── environment.txt (conda list output)
   ```

---

## ✅ **Baseline Validation Script**

Create `scripts/run_baseline.sh`:

```bash
#!/bin/bash
set -e  # Exit on error

echo "=== WASP2 Baseline Pipeline Execution ==="
echo "Started: $(date)"

# Check dependencies
echo "Checking dependencies..."
command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found!"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found!"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found!"; exit 1; }

# Create output directories
mkdir -p baselines/{mapping,counting,analysis}

# Step 1: Counting
echo "Step 1: Counting alleles..."
python -m src.counting count-variants \
    test_data/CD4_ATACseq_Day1_merged_filtered.sort.bam \
    test_data/filter_chr10.vcf \
    --samples NA12878 \
    --region test_data/NA12878_snps_chr10.bed \
    --out baselines/counting/counts.tsv

# Step 2: Analysis
echo "Step 2: Analyzing allelic imbalance..."
python -m src.analysis find-imbalance \
    baselines/counting/counts.tsv \
    --out baselines/analysis/ai_results.tsv

# Save baseline hashes
echo "Saving baseline checksums..."
md5sum baselines/counting/counts.tsv > baselines/counting.md5
md5sum baselines/analysis/ai_results.tsv > baselines/analysis.md5

echo "=== Baseline Complete ==="
echo "Finished: $(date)"

# Display results
echo ""
echo "Results:"
echo "- Counts: $(wc -l < baselines/counting/counts.tsv) rows"
echo "- AI results: $(wc -l < baselines/analysis/ai_results.tsv) rows"
```

---

## 🔄 **Regression Testing Script**

Create `scripts/validate_against_baseline.sh`:

```bash
#!/bin/bash
set -e

echo "=== WASP2 Regression Test ==="

# Run pipeline
./scripts/run_baseline.sh 2>&1 | tee test_run.log

# Compare outputs
echo "Comparing against baseline..."

# Counting comparison
if md5sum -c baselines/counting.md5 2>/dev/null; then
    echo "✓ Counting output matches baseline"
else
    echo "✗ Counting output DIFFERS from baseline!"
    diff baselines/counting/counts.tsv.expected baselines/counting/counts.tsv || true
fi

# Analysis comparison
if md5sum -c baselines/analysis.md5 2>/dev/null; then
    echo "✓ Analysis output matches baseline"
else
    echo "✗ Analysis output DIFFERS from baseline!"
    diff baselines/analysis/ai_results.tsv.expected baselines/analysis/ai_results.tsv || true
fi

echo "=== Regression Test Complete ==="
```

---

## 📋 **Next Steps**

1. **Fix environment.yml** - Add bcftools, bedtools, samtools
2. **Install conda environment**
3. **Run baseline pipeline** (Counting + Analysis minimum)
4. **Save baseline outputs** for regression testing
5. **Document any errors/issues encountered**
6. **(Optional) Run full mapping pipeline** if reference genome available

---

## 🐛 **Known Issues to Watch For**

Based on code review, expect these issues during execution:

1. **Binary search not used** - May be slow on large files
2. **Sample parsing bug** - Check if samples=None causes crash (line 118 __main__.py)
3. **AnnData transpose issue** - If running single-cell variant
4. **Temp file cleanup** - Check if temp files are properly deleted
5. **Error handling** - Look for cryptic error messages (print instead of exceptions)

---

**Document Version**: 1.0
**Last Updated**: 2025-11-15
**Status**: Ready for execution once environment is fixed
