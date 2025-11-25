"""Parity check: run Rust code paths and compare to frozen baselines."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TEST_DATA = ROOT / "test_data"
BAM = TEST_DATA / "CD4_ATACseq_Day1_merged_filtered.sort.bam"
VCF = TEST_DATA / "filter_chr10.vcf"
BED = TEST_DATA / "NA12878_snps_chr10.bed"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

BASELINE_DIR = ROOT / "validation" / "baseline"
CURRENT_RUN_DIR = ROOT / "validation" / "current_run"


def fail(msg: str) -> None:
    print(f"❌ {msg}")
    sys.exit(1)


def run_counts_rust(out_path: Path) -> None:
    env = os.environ.copy()
    env["WASP2_RUST_THREADS"] = env.get("WASP2_RUST_THREADS", "1")
    cmd = [
        "python",
        "-m",
        "src.counting",
        "count-variants",
        str(BAM),
        str(VCF),
        "--regions",
        str(BED),
        "--out_file",
        str(out_path),
    ]
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True, env=env, cwd=ROOT)


def compare_counts(baseline: Path, candidate: Path) -> None:
    b = pd.read_csv(baseline, sep="\t")
    c = pd.read_csv(candidate, sep="\t")
    sort_cols = ["chrom", "pos", "ref", "alt"]
    b = b.sort_values(sort_cols).reset_index(drop=True)
    c = c.sort_values(sort_cols).reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(b, c, check_like=False, check_exact=True)
    except AssertionError as exc:
        fail(f"Counting mismatch vs baseline: {exc}")
    print("✅ Counting matches baseline")


def run_mapping_rust(tmpdir: Path) -> tuple[Path, Path]:
    from validation.generate_mapping_data import write_bams_from_real, collect_qnames
    from mapping.filter_remap_reads import filt_remapped_reads

    tmpdir.mkdir(parents=True, exist_ok=True)
    to_remap, remapped = write_bams_from_real(BAM, tmpdir, n_pairs=50, moved_fraction=0.2)
    out_keep = tmpdir / "keep_rust_t1.bam"
    filt_remapped_reads(str(to_remap), str(remapped), str(out_keep), use_rust=True, threads=1)

    kept_names = collect_qnames(out_keep)
    all_names = collect_qnames(to_remap)
    removed_names = all_names - kept_names

    keep_txt = tmpdir / "mapping_kept_rust.txt"
    removed_txt = tmpdir / "mapping_removed_rust.txt"
    keep_txt.write_text("\n".join(sorted(kept_names)) + "\n")
    removed_txt.write_text("\n".join(sorted(removed_names)) + "\n")
    return keep_txt, removed_txt


def compare_mapping(baseline_keep: Path, baseline_removed: Path, rust_keep: Path, rust_removed: Path) -> None:
    def read_lines(path: Path) -> list[str]:
        return [line.strip() for line in path.read_text().splitlines() if line.strip()]

    b_keep = read_lines(baseline_keep)
    b_removed = read_lines(baseline_removed)
    r_keep = read_lines(rust_keep)
    r_removed = read_lines(rust_removed)

    if b_keep != r_keep:
        fail("Mapping kept read names differ from baseline")
    if b_removed != r_removed:
        fail("Mapping removed read names differ from baseline")
    print("✅ Mapping filter matches baseline (kept/removed names)")


def run_analysis_rust(counts_path: Path, out_path: Path) -> None:
    import wasp2_rust

    results = wasp2_rust.analyze_imbalance(
        str(counts_path),
        min_count=10,
        pseudocount=1,
        method="single",
    )
    df = pd.DataFrame(results)
    df.sort_values(by=df.columns.tolist(), inplace=True)
    df.to_csv(out_path, sep="\t", index=False)


def compare_analysis(baseline: Path, candidate: Path) -> None:
    b = pd.read_csv(baseline, sep="\t")
    c = pd.read_csv(candidate, sep="\t")
    b = b.sort_values(by=b.columns.tolist()).reset_index(drop=True)
    c = c.sort_values(by=c.columns.tolist()).reset_index(drop=True)

    if b.shape != c.shape or list(b.columns) != list(c.columns):
        fail("Analysis output shape/columns differ from baseline")

    numeric_cols = b.select_dtypes(include=["float64", "float32", "int64", "int32", "uint16", "uint32"]).columns
    non_numeric_cols = [col for col in b.columns if col not in numeric_cols]

    if non_numeric_cols and not b[non_numeric_cols].equals(c[non_numeric_cols]):
        fail("Analysis categorical columns differ from baseline")

    for col in numeric_cols:
        if not np.allclose(b[col], c[col], rtol=1e-10, atol=1e-10):
            fail(f"Analysis numeric column mismatch: {col}")

    print("✅ Analysis matches baseline (within tight tolerance)")


def main() -> None:
    shutil.rmtree(CURRENT_RUN_DIR, ignore_errors=True)
    CURRENT_RUN_DIR.mkdir(parents=True, exist_ok=True)

    baseline_counts = BASELINE_DIR / "counts_rust.tsv"
    baseline_analysis = BASELINE_DIR / "analysis_rust.tsv"
    baseline_mapping_keep = BASELINE_DIR / "mapping_kept_rust.txt"
    baseline_mapping_removed = BASELINE_DIR / "mapping_removed_rust.txt"

    for path in [baseline_counts, baseline_analysis, baseline_mapping_keep, baseline_mapping_removed]:
        if not path.exists():
            fail(f"Baseline missing: {path}")

    rust_counts = CURRENT_RUN_DIR / "counts_rust.tsv"
    run_counts_rust(rust_counts)
    compare_counts(baseline_counts, rust_counts)

    mapping_tmp = CURRENT_RUN_DIR / "mapping"
    rust_keep, rust_removed = run_mapping_rust(mapping_tmp)
    compare_mapping(baseline_mapping_keep, baseline_mapping_removed, rust_keep, rust_removed)

    rust_analysis = CURRENT_RUN_DIR / "analysis_rust.tsv"
    run_analysis_rust(rust_counts, rust_analysis)
    compare_analysis(baseline_analysis, rust_analysis)

    shutil.rmtree(mapping_tmp, ignore_errors=True)
    print("\nAll parity checks passed.")


if __name__ == "__main__":
    main()
