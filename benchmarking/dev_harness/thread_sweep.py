#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Iterable, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _ensure_import_paths() -> None:
    root = _repo_root()
    src = root / "src"
    for p in (root, src):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def _parse_int_list(value: str) -> list[int]:
    items: list[int] = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        items.append(int(part))
    return items


def _subset_bam(in_bam: Path, out_bam: Path, regions: list[str]) -> None:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_bam.with_suffix(".tmp.bam")
    cmd = ["samtools", "view", "-b", "-h", str(in_bam), *regions, "-o", str(tmp)]
    subprocess.run(cmd, check=True)
    subprocess.run(["samtools", "sort", "-o", str(out_bam), str(tmp)], check=True)
    subprocess.run(["samtools", "index", str(out_bam)], check=True)
    tmp.unlink(missing_ok=True)


def _subset_vcf(in_vcf: Path, out_vcf_gz: Path, region: str) -> None:
    out_vcf_gz.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["bcftools", "view", "-Oz", "-r", region, "-o", str(out_vcf_gz), str(in_vcf)],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", str(out_vcf_gz)], check=True)


def _count_bam_reads(bam: Path) -> int:
    p = subprocess.run(["samtools", "view", "-c", str(bam)], check=True, capture_output=True, text=True)
    return int(p.stdout.strip())


def _write_bed_from_vcf(
    vcf_gz: Path,
    out_bed: Path,
    sample: str,
    include_indels: bool,
    max_indel_len: int,
) -> None:
    _ensure_import_paths()
    from mapping.intersect_variant_data import vcf_to_bed

    vcf_to_bed(
        vcf_file=vcf_gz,
        out_bed=out_bed,
        samples=[sample],
        include_indels=include_indels,
        max_indel_len=max_indel_len,
    )


def _run_unified(
    subset_bam: Path,
    subset_bed: Path,
    *,
    threads: int,
    bam_threads: int,
    out_dir: Path,
    max_seqs: int,
    indel_mode: bool,
    compress_output: bool,
    compression_threads: int,
) -> dict:
    import wasp2_rust

    out_dir.mkdir(parents=True, exist_ok=True)
    out_r1 = out_dir / "remap_r1.fq"
    out_r2 = out_dir / "remap_r2.fq"

    os.environ["WASP2_BAM_THREADS"] = str(bam_threads)
    os.environ["WASP2_TIMING"] = "1"

    t0 = time.perf_counter()
    stats = wasp2_rust.unified_make_reads_parallel_py(
        str(subset_bam),
        str(subset_bed),
        str(out_r1),
        str(out_r2),
        max_seqs=max_seqs,
        threads=max(1, threads),
        channel_buffer=50_000,
        compression_threads=max(1, compression_threads),
        compress_output=bool(compress_output),
        indel_mode=bool(indel_mode),
        max_indel_size=50,
    )
    wall_s = time.perf_counter() - t0

    # Convert to plain dict (PyO3 returns dict already, but keep explicit)
    result: dict = dict(stats)
    result["wall_s"] = wall_s
    result["threads"] = threads
    result["bam_threads"] = bam_threads
    result["compress_output"] = bool(compress_output)
    result["compression_threads"] = compression_threads
    result["indel_mode"] = bool(indel_mode)
    return result


def main() -> int:
    p = argparse.ArgumentParser(description="Thread sweep for unified pipeline on a BAM/VCF region subset.")
    p.add_argument("--bam", type=Path, required=True, help="Input coordinate-sorted BAM (indexed)")
    p.add_argument("--vcf", type=Path, required=True, help="Input VCF/VCF.GZ (indexed)")
    p.add_argument("--sample", type=str, default="NA12878", help="Sample name for het filtering")
    p.add_argument("--region", type=str, required=True, help="Region (samtools syntax), e.g. chr1:1-2000000")
    p.add_argument("--include-indels", action="store_true", help="Include indels when converting VCF → BED")
    p.add_argument("--max-indel-len", type=int, default=50, help="Max indel length (bp) to include")
    p.add_argument("--threads", type=str, default="1,2,4,8", help="Rayon threads to test (comma-separated)")
    p.add_argument(
        "--bam-threads",
        type=str,
        default="1,2",
        help="htslib BAM decompression threads per worker (comma-separated; WASP2_BAM_THREADS)",
    )
    p.add_argument("--reps", type=int, default=1, help="Replicates per configuration")
    p.add_argument("--max-seqs", type=int, default=64, help="Max haplotype sequences per read pair")
    p.add_argument("--indel-mode", action="store_true", help="Run unified in indel_mode=True (trim combos)")
    p.add_argument("--compress-output", action="store_true", help="Write gzipped FASTQs (slower; for disk storage)")
    p.add_argument("--compression-threads", type=int, default=1, help="gzp threads per FASTQ (if compressing)")
    p.add_argument(
        "--out-tsv",
        type=Path,
        default=Path("benchmarking/dev_harness_out/thread_sweep.tsv"),
        help="Where to write TSV results",
    )
    p.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep the temporary subset + per-run outputs under the output directory",
    )
    args = p.parse_args()

    _ensure_import_paths()

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)

    thread_list = _parse_int_list(args.threads)
    bam_thread_list = _parse_int_list(args.bam_threads)
    if not thread_list or not bam_thread_list:
        raise SystemExit("threads and bam-threads must be non-empty comma-separated lists")

    work_parent: Optional[str] = str(args.out_tsv.parent) if args.keep_tmp else None
    work_dir = Path(tempfile.mkdtemp(prefix="wasp2_thread_sweep_", dir=work_parent))
    try:
        subset_bam = work_dir / "subset.bam"
        subset_vcf = work_dir / "subset.vcf.gz"
        subset_bed = work_dir / "subset.bed"

        regions = [r.strip() for r in args.region.split(",") if r.strip()]
        if not regions:
            raise SystemExit("--region must be non-empty (supports comma-separated regions)")

        print(f"Subsetting BAM to {args.region}...")
        _subset_bam(args.bam, subset_bam, regions)
        subset_reads = _count_bam_reads(subset_bam)
        print(f"Subset BAM reads: {subset_reads:,}")

        print(f"Subsetting VCF to {args.region} and converting to BED (het-only)...")
        _subset_vcf(args.vcf, subset_vcf, args.region)
        _write_bed_from_vcf(
            subset_vcf,
            subset_bed,
            sample=args.sample,
            include_indels=bool(args.include_indels),
            max_indel_len=int(args.max_indel_len),
        )
        bed_variants = sum(1 for _ in subset_bed.open("r", encoding="utf-8"))
        print(f"Subset BED variants: {bed_variants:,}")

        out_fields: list[str] = [
            "threads",
            "bam_threads",
            "rep",
            "wall_s",
            "tree_build_ms",
            "bam_stream_ms",
            "writer_thread_ms",
            "overlap_query_ms",
            "pair_process_ms",
            "send_ms",
            "pairs_processed",
            "pairs_with_variants",
            "pairs_with_snvs_only",
            "pairs_with_indels_only",
            "pairs_with_snvs_and_indels",
            "haplotypes_written",
            "pairs_kept",
            "pairs_keep_no_flip",
            "pairs_skipped_unmappable",
            "pairs_haplotype_failed",
            "orphan_reads",
            "compress_output",
            "compression_threads",
            "indel_mode",
        ]

        with args.out_tsv.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
            w.writeheader()
            for rep in range(args.reps):
                for t in thread_list:
                    for bt in bam_thread_list:
                        run_dir = work_dir / f"run_t{t}_bt{bt}_rep{rep}"
                        res = _run_unified(
                            subset_bam,
                            subset_bed,
                            threads=t,
                            bam_threads=bt,
                            out_dir=run_dir,
                            max_seqs=int(args.max_seqs),
                            indel_mode=bool(args.indel_mode),
                            compress_output=bool(args.compress_output),
                            compression_threads=int(args.compression_threads),
                        )
                        row = {k: res.get(k) for k in out_fields if k not in {"rep"}}
                        row["rep"] = rep
                        w.writerow(row)

                        # Keep disk usage manageable unless explicitly debugging.
                        if not args.keep_tmp:
                            shutil.rmtree(run_dir, ignore_errors=True)

        print(f"Wrote: {args.out_tsv}")
        return 0
    finally:
        if args.keep_tmp:
            print(f"Keeping work dir: {work_dir}")
        else:
            shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
