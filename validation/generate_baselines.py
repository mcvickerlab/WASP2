"""Generate frozen baselines for parity checks (Rust paths)."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = ROOT / "test_data"
BAM = TEST_DATA / "CD4_ATACseq_Day1_merged_filtered.sort.bam"
VCF = TEST_DATA / "filter_chr10.vcf"
BED = TEST_DATA / "NA12878_snps_chr10.bed"

BASELINE_DIR = ROOT / "validation" / "baseline"
CURRENT_RUN_DIR = ROOT / "validation" / "current_run"

# Ensure local src on path for module imports
import sys

if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))


def ensure_dirs() -> None:
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    CURRENT_RUN_DIR.mkdir(parents=True, exist_ok=True)


def run_counts(out_path: Path) -> None:
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
    subprocess.run(cmd, check=True, env=env, cwd=ROOT)


def run_mapping(tmpdir: Path) -> tuple[Path, Path]:
    """Generate kept/removed names using the Rust mapping filter."""
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


def run_analysis(counts_path: Path, out_path: Path) -> None:
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


def main() -> None:
    ensure_dirs()
    shutil.rmtree(CURRENT_RUN_DIR, ignore_errors=True)
    CURRENT_RUN_DIR.mkdir(parents=True, exist_ok=True)

    counts_path = CURRENT_RUN_DIR / "counts_rust.tsv"
    run_counts(counts_path)
    shutil.copy(counts_path, BASELINE_DIR / "counts_rust.tsv")

    mapping_tmp = CURRENT_RUN_DIR / "mapping"
    rust_keep, rust_removed = run_mapping(mapping_tmp)
    shutil.copy(rust_keep, BASELINE_DIR / "mapping_kept_rust.txt")
    shutil.copy(rust_removed, BASELINE_DIR / "mapping_removed_rust.txt")

    analysis_path = CURRENT_RUN_DIR / "analysis_rust.tsv"
    run_analysis(counts_path, analysis_path)
    shutil.copy(analysis_path, BASELINE_DIR / "analysis_rust.tsv")

    shutil.rmtree(mapping_tmp, ignore_errors=True)
    print("Baselines written to", BASELINE_DIR)


if __name__ == "__main__":
    main()
