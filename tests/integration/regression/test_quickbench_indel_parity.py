import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"

for p in (ROOT, SRC):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


@pytest.mark.unit
def test_quickbench_indel_parity(tmp_path: Path) -> None:
    """Unified make-reads matches the multi-pass path on a simple INDEL dataset (no trim combos)."""
    wasp2_rust = pytest.importorskip("wasp2_rust")

    from benchmarking.quickbench.fastq_utils import counter_diff, fastq_counter
    from benchmarking.quickbench.synthetic_dataset import (
        quickbench_indel_variants,
        write_bed,
        write_synthetic_bam_indel,
    )
    from mapping.intersect_variant_data import intersect_reads, process_bam
    from mapping.make_remap_reads import write_remap_bam

    bam = tmp_path / "synthetic_indel.bam"
    bed = tmp_path / "variants_indel.bed"
    write_synthetic_bam_indel(bam)
    write_bed(bed, quickbench_indel_variants())

    baseline_dir = tmp_path / "baseline"
    baseline_dir.mkdir()
    to_remap_bam = baseline_dir / "to_remap.bam"
    keep_bam = baseline_dir / "keep.bam"
    remap_reads_txt = baseline_dir / "remap_reads.txt"
    intersect_bed = baseline_dir / "intersect.bed"
    baseline_r1 = baseline_dir / "baseline_r1.fq"
    baseline_r2 = baseline_dir / "baseline_r2.fq"

    process_bam(
        bam_file=str(bam),
        vcf_bed=str(bed),
        remap_bam=str(to_remap_bam),
        remap_reads=str(remap_reads_txt),
        keep_bam=str(keep_bam),
        is_paired=True,
        threads=1,
    )
    intersect_reads(
        remap_bam=str(to_remap_bam),
        vcf_bed=str(bed),
        out_bed=str(intersect_bed),
        num_samples=1,
    )
    write_remap_bam(
        bam_file=str(to_remap_bam),
        intersect_file=str(intersect_bed),
        r1_out=str(baseline_r1),
        r2_out=str(baseline_r2),
        samples=["SYNTH"],
        max_seqs=64,
        include_indels=True,
    )

    unified_dir = tmp_path / "unified"
    unified_dir.mkdir()
    unified_r1 = unified_dir / "unified_r1.fq"
    unified_r2 = unified_dir / "unified_r2.fq"

    wasp2_rust.unified_make_reads_py(
        str(bam),
        str(bed),
        str(unified_r1),
        str(unified_r2),
        max_seqs=64,
        threads=1,
        compression_threads=1,
        compress_output=False,
        indel_mode=False,
    )

    baseline_counter = fastq_counter(baseline_r1, baseline_r2)
    unified_counter = fastq_counter(unified_r1, unified_r2)
    only_baseline, only_unified = counter_diff(baseline_counter, unified_counter)

    assert only_baseline == [] and only_unified == [], (
        "INDEL parity mismatch between multi-pass and unified outputs.\n"
        f"Only in baseline: {only_baseline[:5]}\n"
        f"Only in unified:  {only_unified[:5]}"
    )

