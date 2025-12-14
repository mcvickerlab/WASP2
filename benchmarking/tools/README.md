# Benchmarking tools

Small helper scripts for reproducible profiling and parameter sweeps.

## `log_unified_stats.py`

Prints a single greppable line from a `unified_stats.json` produced by the Rust unified pipeline.

Example:

```bash
python benchmarking/tools/log_unified_stats.py /path/to/unified_stats.json --label atac_snv_step2
```

## `make_hg00731_subset.sh`

Creates deterministic HG00731 RNA-seq subsets (chr15–22 + 2%/10% subsamples) for quick profiling runs.

```bash
benchmarking/tools/make_hg00731_subset.sh --out-dir /tmp/hg00731_subsets
```

## `thread_sweep_unified_make_reads.sh`

Runs a thread + compression sweep for `unified_make_reads_parallel_py` on an existing BAM+BED.

```bash
benchmarking/tools/thread_sweep_unified_make_reads.sh \
  --bam /tmp/hg00731_subsets/HG00731_chr15-22_sub10pct.bam \
  --bed /tmp/hg00731_subsets/HG00731_chr15-22_het_only.bed \
  --out-dir /tmp/hg00731_thread_sweep
```
