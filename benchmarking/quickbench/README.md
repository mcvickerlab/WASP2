# Quickbench

Fast, deterministic checks used to keep WASP2-Rust optimizations honest.

This is **not** the 150M-read benchmark suite. It’s designed to run in seconds.

## SNV Parity (Unified vs Multi-pass)

Runs a tiny synthetic paired-end BAM through:

1. **Multi-pass** pipeline (`process_bam` → `intersect_reads` → `remap_all_chromosomes`)
2. **Unified** pipeline (`wasp2_rust.unified_make_reads_py`)

Then compares the emitted remap FASTQs as an **order-independent multiset**
of `(orig_name, mate, sequence, qualities)`.

Run:

```bash
python benchmarking/quickbench/run_quickbench.py snv-parity
```

Outputs:

- `benchmarking/quickbench_out/quickbench_snv_parity.json`

To keep intermediate files for debugging:

```bash
python benchmarking/quickbench/run_quickbench.py snv-parity --keep-tmp
```

