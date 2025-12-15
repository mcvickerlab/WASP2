# Dev Harness

Tools for fast iteration on real data without running the full 150M benchmarks.

## Compare Unified vs Multi-pass on a Region

This subsets both the BAM and BED to a region, then runs:

- **Unified** (`unified_make_reads_parallel_py`, `indel_mode=False`)
- **Multi-pass** (`process_bam` → `intersect_reads` → `remap_all_chromosomes`)

and compares the emitted remap FASTQs as an order-independent multiset of
`(orig_name, mate, sequence, qualities)`.

Run:

```bash
python benchmarking/dev_harness/compare_subset.py \
  --bam /path/to/input.bam \
  --bed /path/to/variants.bed \
  --region chr1:1,000,000-2,000,000
```

Outputs:
- `benchmarking/dev_harness_out/subset_compare_summary.json`

Notes:
- `--indel-mode` runs unified trim-combo mode; parity compare is disabled because the multi-pass baseline does not generate trim combos.

## Thread Sweep (Unified Only)

This runs unified on a subset BAM/BED repeatedly with different:
- Rayon worker threads (`--threads`)
- htslib BAM decompression threads per worker (`--bam-threads`, via `WASP2_BAM_THREADS`)

Run:

```bash
python benchmarking/dev_harness/thread_sweep.py \
  --bam /path/to/input.bam \
  --vcf /path/to/variants.vcf.gz \
  --sample NA12878 \
  --region chr1:1-2000000 \
  --threads 1,2,4,8,16 \
  --bam-threads 1,2,4 \
  --out-tsv benchmarking/dev_harness_out/thread_sweep.tsv
```
