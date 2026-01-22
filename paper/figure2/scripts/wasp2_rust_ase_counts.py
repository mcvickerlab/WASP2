#!/usr/bin/env python3
"""
Fast per-variant ref/alt counting for Figure 2 using `wasp2_rust.BamCounter`.

Outputs TSV:
  chrom  pos  ref  alt  ref_count  alt_count  other_count

Only biallelic heterozygous SNPs (len(ref)==1 and len(alt)==1) are counted.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pysam
from wasp2_rust import BamCounter


def iter_het_snps(vcf_path: Path, sample: str):
    vcf = pysam.VariantFile(str(vcf_path))
    if sample not in vcf.header.samples:
        raise SystemExit(
            f"Sample '{sample}' not found in VCF header samples: {list(vcf.header.samples)}"
        )

    for rec in vcf:
        if not rec.alts or len(rec.alts) != 1:
            continue
        ref = rec.ref
        alt = rec.alts[0]
        if len(ref) != 1 or len(alt) != 1:
            continue
        gt = rec.samples[sample].get("GT")
        if gt not in {(0, 1), (1, 0)}:
            continue
        yield rec.chrom, int(rec.pos), ref, alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--min-baseq", type=int, default=0)
    ap.add_argument(
        "--max-sites", type=int, default=0, help="Optional cap for quick tests (0 = no cap)"
    )
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    regions = []
    for i, (chrom, pos, ref, alt) in enumerate(iter_het_snps(args.vcf, args.sample), start=1):
        regions.append((chrom, pos, ref, alt))
        if args.max_sites and i >= args.max_sites:
            break

    counter = BamCounter(str(args.bam))
    counts = counter.count_alleles(
        regions,
        min_qual=int(args.min_baseq),
        min_mapq=int(args.min_mapq),
        threads=max(1, int(args.threads)),
    )

    with args.out.open("w") as f:
        f.write("chrom\tpos\tref\talt\tref_count\talt_count\tother_count\n")
        for (chrom, pos, ref, alt), (ref_c, alt_c, other_c) in zip(regions, counts):
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{ref_c}\t{alt_c}\t{other_c}\n")


if __name__ == "__main__":
    main()
