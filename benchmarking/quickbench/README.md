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

## INDEL Parity (No Trim)

This validates that **allele substitution** (including indels) matches between:
- multi-pass remap (`process_bam` → `intersect_reads` → `remap_all_chromosomes`)
- unified remap (`unified_make_reads_py`) with `indel_mode=False`

Run:

```bash
python benchmarking/quickbench/run_quickbench.py indel-parity
```

## INDEL Trim Invariants (PI-style combos)

This validates the **N+1 trim-combo behavior** for a +2bp insertion, ensuring:
- output reads stay at the original read length
- exactly 3 trimmed alternates are emitted for the insertion haplotype
- `_WASP_..._total_seqs` matches the actual emitted count

Run:

```bash
python benchmarking/quickbench/run_quickbench.py indel-trim-invariants
```
