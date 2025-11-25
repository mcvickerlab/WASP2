# Validation harness for the Rust rewrite

Goal: keep a frozen, reproducible baseline from the current pipeline and provide a one-command parity check as we remove Python fallbacks.

## What’s included
- `generate_baselines.py` – builds deterministic baseline outputs using the current code path (Rust by default):
  - Counting: `validation/baseline/counts_rust.tsv`
  - Analysis: `validation/baseline/analysis_rust.tsv`
  - Mapping filter: `validation/baseline/mapping_kept_rust.txt` and `validation/baseline/mapping_removed_rust.txt`
- `compare_to_baseline.py` – runs the Rust path and diffs against the baselines. Fails (exit 1) on any mismatch.

## Usage
1) Activate the WASP2 environment and ensure the Rust extension is built:
   ```bash
   conda activate WASP2
   export LIBCLANG_PATH=$CONDA_PREFIX/lib
   export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
   export BINDGEN_EXTRA_CLANG_ARGS="-I/usr/include"
   (cd rust && maturin develop --release)
   ```
2) Generate baselines once (creates `validation/baseline/*`):
   ```bash
   python validation/generate_baselines.py
   ```
3) Run parity check after changes:
   ```bash
   python validation/compare_to_baseline.py
   ```
   - Counting parity: exact TSV match.
   - Mapping filter parity: identical kept/removed read names.
   - Analysis parity: floating-point compare with tight tolerance (1e-10 relative/absolute).

## Notes
- Baselines use the small chr10 test bundle in `test_data/` and a deterministic synthetic remap set (20k read pairs, 20% moved).
- Temporary files live under `validation/current_run/` and are ignored by git.
