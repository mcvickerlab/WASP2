"""
Rust vs Python parity test for allele counting.

Reimplements the Rust BamCounter algorithm in pure Python + pysam and
compares allele counts at every variant position. The Python reference
matches the Rust semantics exactly:

- Single BAM fetch per chromosome spanning all variant positions
- Each read assigned to the earliest-encounter-index SNP it overlaps
- Seen-reads set accumulates across positions (each read counted once)
- BAM flag filtering: unmapped, secondary, supplementary, duplicate skipped

This catches any numerical differences between the Rust and Python
implementations that golden-file tests (which compare Rust output to
Rust-generated baselines) would miss.
"""

from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[1]
SHARED_DATA = ROOT / "tests" / "shared_data"
SANITY_DATA = ROOT / "tests" / "sanity" / "data"


# ---------------------------------------------------------------------------
# Python counting algorithm matching Rust BamCounter semantics
# ---------------------------------------------------------------------------
def python_count_snp_alleles(bam_path: str, chrom: str, snp_list: list[tuple]):
    """Count ref/alt/other alleles using pure Python + pysam.

    Reimplements the Rust BamCounter algorithm (bam_counter.rs:228-370) in
    Python to enable exact parity comparison. Key semantics matched:

    1. Single fetch for the entire chromosome span (min_pos to max_pos).
    2. Each read is assigned to the earliest-encounter-index SNP it has an
       aligned base at. A read is counted at exactly one position.
    3. Seen-reads set accumulates across all positions (paired-end mates
       and multi-overlap reads are counted only once per chromosome).
    4. BAM flag filtering: skip unmapped, secondary, supplementary, duplicate.

    Args:
        bam_path: Path to BAM file.
        chrom: Chromosome name.
        snp_list: List of (pos_1based, ref, alt) tuples, in encounter order.

    Returns:
        List of (chrom, pos_1based, ref_count, alt_count, other_count) tuples.
    """
    if not snp_list:
        return []

    # Build position -> list of (encounter_index, ref, alt)
    pos_map: dict[int, list[tuple]] = {}
    for enc_idx, (pos, ref, alt) in enumerate(snp_list):
        pos_map.setdefault(pos, []).append((enc_idx, ref, alt))

    # Initialize counts per encounter index
    counts: dict[int, list[int]] = {i: [0, 0, 0] for i in range(len(snp_list))}

    min_pos = min(pos for pos, _, _ in snp_list)
    max_pos = max(pos for pos, _, _ in snp_list)

    seen_reads: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, max(0, min_pos - 1), max_pos + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue

            qname = read.query_name
            if qname in seen_reads:
                continue

            seq = read.query_sequence
            if seq is None:
                continue

            # Find earliest-encounter-index SNP this read has an aligned base at
            best = None  # (encounter_idx, ref, alt, qpos)
            for qpos, refpos in read.get_aligned_pairs(True):
                if qpos is None or refpos is None:
                    continue
                pos1 = refpos + 1
                if pos1 in pos_map:
                    for enc_idx, ref, alt in pos_map[pos1]:
                        if best is None or enc_idx < best[0]:
                            best = (enc_idx, ref, alt, qpos)

            if best is not None:
                enc_idx, ref, alt, qpos = best
                base = seq[qpos]
                if base == ref:
                    counts[enc_idx][0] += 1
                elif base == alt:
                    counts[enc_idx][1] += 1
                else:
                    counts[enc_idx][2] += 1
                seen_reads.add(qname)

    return [
        (chrom, snp_list[i][0], counts[i][0], counts[i][1], counts[i][2])
        for i in range(len(snp_list))
    ]


def parse_het_variants_from_vcf(vcf_path: str, sample: str | None = None):
    """Extract heterozygous variant positions from a VCF file.

    Returns:
        Dict mapping chrom -> list of (pos_1based, ref, alt).
    """
    variants_by_chrom: dict[str, list[tuple]] = {}

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            # Get genotype for first sample or specified sample
            if sample and sample in rec.samples:
                gt = rec.samples[sample]["GT"]
            else:
                gt = rec.samples[list(rec.samples)[0]]["GT"]

            if gt == (0, 1) or gt == (1, 0):
                chrom = rec.chrom
                if chrom not in variants_by_chrom:
                    variants_by_chrom[chrom] = []
                variants_by_chrom[chrom].append((rec.pos, rec.ref, rec.alts[0]))

    return variants_by_chrom


# ---------------------------------------------------------------------------
# Parity tests
# ---------------------------------------------------------------------------
class TestCountingParity:
    """Compare Python and Rust counting on the same data."""

    @pytest.fixture
    def shared_data(self):
        if not SHARED_DATA.exists() or not (SHARED_DATA / "sample1.bam").exists():
            pytest.skip("Shared test data not available")
        return {
            "bam": SHARED_DATA / "sample1.bam",
            "vcf": SHARED_DATA / "variants.vcf.gz",
        }

    def test_counting_parity_shared_data(self, shared_data):
        """Rust and Python must produce identical allele counts on shared test data."""
        import wasp2_rust

        bam_path = str(shared_data["bam"])
        vcf_path = str(shared_data["vcf"])

        # Get het variants by chromosome
        variants_by_chrom = parse_het_variants_from_vcf(vcf_path, sample="sample1")

        # Collect Python counts
        python_counts = {}
        for chrom, snp_list in variants_by_chrom.items():
            results = python_count_snp_alleles(bam_path, chrom, snp_list)
            for chrom_r, pos, ref_c, alt_c, other_c in results:
                python_counts[(chrom_r, pos)] = (ref_c, alt_c, other_c)

        # Collect Rust counts
        all_regions = []
        for chrom, snp_list in variants_by_chrom.items():
            for pos, ref, alt in snp_list:
                all_regions.append((chrom, pos, ref, alt))

        counter = wasp2_rust.BamCounter(bam_path)
        rust_results = counter.count_alleles(all_regions, min_qual=0, threads=1)

        rust_counts = {}
        for (chrom, pos, _ref, _alt), (ref_c, alt_c, other_c) in zip(all_regions, rust_results):
            rust_counts[(chrom, pos)] = (ref_c, alt_c, other_c)

        # Compare
        assert len(python_counts) > 0, "Should have variants to compare"
        assert len(python_counts) == len(rust_counts), (
            f"Variant count mismatch: Python={len(python_counts)}, Rust={len(rust_counts)}"
        )

        mismatches = []
        for key in sorted(python_counts):
            py = python_counts[key]
            rs = rust_counts.get(key)
            if rs is None:
                mismatches.append(f"{key}: missing from Rust results")
            elif py != rs:
                mismatches.append(
                    f"{key[0]}:{key[1]} - Python=({py[0]},{py[1]},{py[2]}) "
                    f"Rust=({rs[0]},{rs[1]},{rs[2]})"
                )

        assert not mismatches, (
            f"Counting parity failures ({len(mismatches)}/{len(python_counts)}):\n"
            + "\n".join(mismatches[:20])
        )


class TestCountingParitySanity:
    """Compare Python and Rust counting on chr21 real data (larger dataset)."""

    @pytest.fixture
    def sanity_data(self):
        if not SANITY_DATA.exists() or not (SANITY_DATA / "chr21.bam").exists():
            pytest.skip("Sanity test data not available (run 'make download-sanity-data')")
        return {
            "bam": SANITY_DATA / "chr21.bam",
            "vcf": SANITY_DATA / "chr21.vcf.gz",
        }

    @pytest.mark.sanity
    def test_counting_parity_chr21(self, sanity_data):
        """Rust and Python must produce identical counts on real chr21 data.

        This is the definitive parity test — 33k+ variants, real sequencing
        data from HG00731. Any systematic bias would surface here.
        """
        import wasp2_rust

        bam_path = str(sanity_data["bam"])
        vcf_path = str(sanity_data["vcf"])

        variants_by_chrom = parse_het_variants_from_vcf(vcf_path)

        # Python counts
        python_counts = {}
        for chrom, snp_list in variants_by_chrom.items():
            results = python_count_snp_alleles(bam_path, chrom, snp_list)
            for chrom_r, pos, ref_c, alt_c, other_c in results:
                python_counts[(chrom_r, pos)] = (ref_c, alt_c, other_c)

        # Rust counts
        all_regions = []
        for chrom, snp_list in variants_by_chrom.items():
            for pos, ref, alt in snp_list:
                all_regions.append((chrom, pos, ref, alt))

        counter = wasp2_rust.BamCounter(bam_path)
        rust_results = counter.count_alleles(all_regions, min_qual=0, threads=1)

        rust_counts = {}
        for (chrom, pos, _ref, _alt), (ref_c, alt_c, other_c) in zip(all_regions, rust_results):
            rust_counts[(chrom, pos)] = (ref_c, alt_c, other_c)

        # Compare
        total = len(python_counts)
        assert total > 1000, f"Expected 1000+ variants, got {total}"

        mismatches = []
        for key in sorted(python_counts):
            py = python_counts[key]
            rs = rust_counts.get(key)
            if rs is None:
                mismatches.append(f"{key}: missing from Rust")
            elif py != rs:
                mismatches.append(
                    f"{key[0]}:{key[1]} - Python=({py[0]},{py[1]},{py[2]}) "
                    f"Rust=({rs[0]},{rs[1]},{rs[2]})"
                )

        mismatch_rate = len(mismatches) / total if total > 0 else 0
        assert not mismatches, (
            f"Counting parity failures: {len(mismatches)}/{total} "
            f"({mismatch_rate:.1%}):\n"
            + "\n".join(mismatches[:20])
            + (f"\n... and {len(mismatches) - 20} more" if len(mismatches) > 20 else "")
        )
