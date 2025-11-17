# WASP2 Test Data Bundle

This directory contains NA12878 chr10 test data for regression testing.

## Files:
- `filter_chr10.vcf`: NA12878 chr10 SNPs
- `NA12878_snps_chr10.bed`: matching intervals
- `CD4_ATACseq_Day1_merged_filtered.sort.bam` + `.bai`: small BAM subset
- `as_counts.txt`: toy ASE counts
- `genome.fa` + BWA indices: chr10 reference genome
- `wasp2_test_bundle.tar.gz`: compressed archive of all test data

These files are used by the automated regression test suite in `tests/regression/`.
