# WASP2 Nature Methods Benchmark Summary

## Dataset
- **Sample**: HG00731 (1000 Genomes, Puerto Rican)
- **Data**: ERR1050079 RNA-seq (56M paired-end reads)
- **Variants**: 2,206,785 heterozygous SNPs from 1000 Genomes VCF
- **Threads**: 8

---

## Speed Benchmark Results

### WASP2 Pipeline (Rust)

| Step | Time (s) | Description |
|------|----------|-------------|
| STAR initial alignment | 475.9 | Align 56M reads to GRCh38 |
| **WASP2 unified make-reads** | **16.9** | Identify reads overlapping variants, generate swapped reads |
| STAR remap | 48.0 | Realign swapped reads |
| **WASP2 filter** | **13.0** | Filter remapped reads |
| **WASP2 total** | **29.8** | Just WASP steps |
| **Full pipeline** | **555.8** | Including STAR |

### WASP2 Allele Counting (Rust BamCounter)

| Metric | Value |
|--------|-------|
| Variants processed | 2,206,785 |
| Time | 15.73s |
| **Speed** | **140,329 variants/sec** |

### Comparison to Published Methods

| Method | WASP Time | Source |
|--------|-----------|--------|
| WASP1 (Python) | 1541s (25:41) | STAR+WASP paper |
| STAR+WASP builtin | 332s (5:32) | STAR+WASP paper |
| **WASP2-Rust** | **30s** | This benchmark |

**Speedup vs WASP1: 51x**
**Speedup vs STAR+WASP: 11x** (WASP steps only)

---

## Accuracy Results

### WASP2 SNP Counting

| Metric | Value |
|--------|-------|
| Variants with coverage | 156,457 / 2,206,785 (7.1%) |
| Total reads counted | 1,675,780 |
| Mean coverage | 10.7 reads/variant |
| Mean REF ratio | 0.520 |

**Mean REF ratio of 0.520 indicates minimal reference bias** (expected 0.50 for unbiased).

---

## Key Findings

1. **WASP2-Rust is 51x faster than WASP1** for the same filtering task
2. **WASP2-Rust is 11x faster than STAR+WASP builtin** (WASP steps only)
3. **WASP2 allele counting achieves 140K variants/second**
4. **Reference bias is minimal** (mean REF ratio = 0.52)

---

## Competitor Tools

### GATK ASEReadCounter
- **Limitation**: Outputs per-base, not per-variant (problematic for INDELs)
- **INDEL support**: Cannot properly count INDEL alleles

### phASER
- **Focus**: Haplotype-level ASE with read-backed phasing
- **INDEL support**: Limited

### biastools
- **Focus**: Reference bias quantification
- **INDEL support**: Yes, via CIGAR analysis
- **Speed**: Python-based, slower than Rust

---

## Files Generated

```
benchmarking/nature_methods_results/
├── wasp2_snp_counts.csv     # 2.2M SNP allele counts
├── timing.csv               # Speed measurements
└── BENCHMARK_SUMMARY.md     # This file

benchmarking/star_wasp_comparison/results/unified_*/
├── benchmark_results.json   # Pipeline timing
├── A_sorted.bam             # Aligned reads
└── remap_keep.bam           # WASP-filtered reads
```
