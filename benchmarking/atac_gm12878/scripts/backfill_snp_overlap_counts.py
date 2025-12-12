#!/usr/bin/env python3
"""
Backfill baseline-filtered variant-overlap denominators for SNP-only ATAC runs.

Panel C uses "Original" counts defined as reads entering the WASP remap+filter
path after baseline BAM filters. INDEL runs already emit split overlap keys.
SNP-only runs do not, so we recompute the denominator using the unified pipeline
on the saved to_remap BAMs and write split keys into their benchmark JSONs.

This script is safe to rerun: it overwrites only the split overlap keys.
"""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]

PY_RUN_DIR = REPO_ROOT / "benchmarking/atac_gm12878/results/wasp2python_snp_fixed_2025-12-07_13-12-00"
RUST_RUN_DIR = REPO_ROOT / "benchmarking/atac_gm12878/results/wasp2rust_snp_fixed_2025-12-10_00-45-32"

AHO_VCF = Path("/iblm/netapp/data1/aho/variants/NA12878.vcf.gz")


def build_bed(vcf_path: Path, out_bed: Path, sample: str | None) -> None:
    import wasp2_rust

    samples = [sample] if sample else None
    wasp2_rust.vcf_to_bed_py(
        str(vcf_path),
        str(out_bed),
        samples=samples,
        het_only=True,
        include_indels=False,
        max_indel_len=0,
        include_genotypes=True,
    )


def unified_pairs_with_variants(to_remap_bam: Path, bed: Path, tmp_dir: Path) -> int:
    from wasp2_rust import unified_make_reads_parallel_py

    r1 = tmp_dir / "remap_r1.fq"
    r2 = tmp_dir / "remap_r2.fq"
    stats = unified_make_reads_parallel_py(
        str(to_remap_bam),
        str(bed),
        str(r1),
        str(r2),
        max_seqs=64,
        threads=8,
        compress_output=False,
        indel_mode=False,
    )
    return int(stats.get("pairs_with_variants", 0))


def patch_json(run_dir: Path, pairs_with_variants: int) -> None:
    json_path = run_dir / "benchmark_results.json"
    d = json.loads(json_path.read_text())

    pre_reads = pairs_with_variants * 2
    post_reads = int(d.get("remap_keep_reads", 0))

    d.update(
        {
            "snv_only_overlap_reads_pre": pre_reads,
            "indel_only_overlap_reads_pre": 0,
            "snv_indel_overlap_reads_pre": 0,
            "snv_only_overlap_reads_post": post_reads,
            "indel_only_overlap_reads_post": 0,
            "snv_indel_overlap_reads_post": 0,
            "pre_definition": "baseline_filtered_variant_overlaps",
        }
    )
    json_path.write_text(json.dumps(d, indent=2) + "\n")

    old_pre = int(d.get("original_reads", 0)) - int(d.get("keep_reads", 0))
    print(f"{run_dir.name}: old pre={old_pre:,} reads, new pre={pre_reads:,} reads, post={post_reads:,}")


def backfill_one(run_dir: Path, bam_name: str, vcf_path: Path, sample: str | None) -> None:
    to_remap_bam = run_dir / bam_name
    if not to_remap_bam.exists():
        raise FileNotFoundError(f"Missing to_remap BAM: {to_remap_bam}")

    with tempfile.TemporaryDirectory(prefix="wasp2_backfill_", dir=str(run_dir)) as td:
        tmp_dir = Path(td)
        bed_path = tmp_dir / "variants_snp_only.bed"
        build_bed(vcf_path, bed_path, sample=sample)
        pairs = unified_pairs_with_variants(to_remap_bam, bed_path, tmp_dir)
        patch_json(run_dir, pairs)

        # Cleanup any outputs we created besides temp dir (should all be inside tmp)
        for p in tmp_dir.glob("*"):
            try:
                if p.is_file():
                    p.unlink()
                else:
                    shutil.rmtree(p, ignore_errors=True)
            except Exception:
                pass


def main() -> None:
    print("Backfilling SNP-only ATAC runs with baseline-filtered denominators...")

    # Python SNP run uses a local het-only VCF already filtered for NA12878.
    py_vcf = PY_RUN_DIR / "het_only.vcf.gz"
    if not py_vcf.exists():
        py_vcf = AHO_VCF
        py_sample = "NA12878"
    else:
        py_sample = None
    backfill_one(PY_RUN_DIR, "original_to_remap.bam", py_vcf, py_sample)

    # Rust SNP run: use Aho NA12878 VCF with sample filtering.
    backfill_one(RUST_RUN_DIR, "to_remap_actual.bam", AHO_VCF, "NA12878")

    print("Done.")


if __name__ == "__main__":
    main()

