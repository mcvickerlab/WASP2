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

