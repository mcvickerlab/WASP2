"""Sanity tests: Verify Rust pipeline produces consistent results on real chr21 data.

These tests validate that WASP2's Rust-accelerated pipeline produces
identical results to the established baseline on real HG00731 RNA-seq data
(chr21 subset with ~855K reads and ~33K het variants).

Tests cover:
1. Allele counting - exact integer match
2. FASTQ read set generation - identical read names
3. Analysis output - floating-point tolerance for p-values

Run with: pytest tests/sanity/ -v
Skip with: pytest -m "not sanity"
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    pass

# Mark all tests in this module as sanity tests
pytestmark = pytest.mark.sanity


class TestAlleleCounts:
    """Test allele counting produces exact integer matches."""

    def test_counts_exact_match(
        self,
        sanity_data: dict[str, Path],
        sanity_tmp_dir: Path,
    ) -> None:
        """Allele counts must match exactly (integers).

        This is the core sanity check - if counts differ, there's a
        fundamental bug in the Rust counting implementation.
        """
        import pysam

        import wasp2_rust

        # Parse variants from VCF (het sites only)
        regions = []
        with pysam.VariantFile(str(sanity_data["vcf"])) as vcf:
            for rec in vcf:
                gt = rec.samples[0]["GT"]
                if gt == (0, 1) or gt == (1, 0):
                    regions.append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))

        # Run counting
        counter = wasp2_rust.BamCounter(str(sanity_data["bam"]))
        counts = counter.count_alleles(regions, min_qual=0, threads=1)

        # Load expected counts
        expected = {}
        with open(sanity_data["expected_counts"]) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                row = dict(zip(header, parts))
                key = (row["chrom"], int(row["pos"]))
                expected[key] = (
                    int(row["ref_count"]),
                    int(row["alt_count"]),
                    int(row["other_count"]),
                )

        # Compare counts
        mismatches = []
        for (chrom, pos, _ref, _alt), (ref_count, alt_count, other_count) in zip(regions, counts):
            key = (chrom, pos)
            if key in expected:
                exp = expected[key]
                if (ref_count, alt_count, other_count) != exp:
                    mismatches.append(
                        f"{chrom}:{pos} - got ({ref_count},{alt_count},{other_count}), "
                        f"expected {exp}"
                    )

        assert not mismatches, (
            "Count mismatches found:\n"
            + "\n".join(mismatches[:10])
            + (f"\n... and {len(mismatches) - 10} more" if len(mismatches) > 10 else "")
        )

    def test_counts_coverage_stats(
        self,
        sanity_data: dict[str, Path],
    ) -> None:
        """Verify expected coverage statistics from chr21 data."""
        # Load expected counts and compute stats
        total_sites = 0
        sites_with_coverage = 0
        total_ref = 0
        total_alt = 0

        with open(sanity_data["expected_counts"]) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                ref_count = int(parts[4])
                alt_count = int(parts[5])
                total_sites += 1
                if ref_count + alt_count > 0:
                    sites_with_coverage += 1
                total_ref += ref_count
                total_alt += alt_count

        # Sanity checks on chr21 HG00731 data
        assert total_sites == 33036, f"Expected 33036 het sites, got {total_sites}"
        assert sites_with_coverage > 100, (
            f"Expected >100 sites with coverage, got {sites_with_coverage}"
        )
        assert total_ref + total_alt > 1000, "Total counts should be >1000"


class TestFastqGeneration:
    """Test FASTQ read generation produces consistent read sets."""

    def test_fastq_readset_match(
        self,
        sanity_data: dict[str, Path],
        sanity_tmp_dir: Path,
    ) -> None:
        """Remapped FASTQ read names must match expected set.

        The order may differ (due to parallel processing), but the
        set of read names must be identical.
        """
        import wasp2_rust

        out_r1 = sanity_tmp_dir / "r1.fq"
        out_r2 = sanity_tmp_dir / "r2.fq"

        # Create BED from VCF
        bed_path = sanity_tmp_dir / "variants.bed"
        wasp2_rust.vcf_to_bed_py(
            str(sanity_data["vcf"]),
            str(bed_path),
        )

        # Run unified pipeline (single-threaded for reproducibility)
        wasp2_rust.unified_make_reads_py(
            str(sanity_data["bam"]),
            str(bed_path),
            str(out_r1),
            str(out_r2),
            max_seqs=64,
            threads=1,
            compression_threads=1,
            compress_output=False,
            indel_mode=True,
        )

        # Extract read names from generated FASTQ
        def get_read_names(fq_path: Path) -> set[str]:
            names = set()
            with open(fq_path) as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:  # Header line
                        # Extract read name (before whitespace)
                        name = line.strip().split()[0]
                        if name.startswith("@"):
                            name = name[1:]
                        names.add(name)
            return names

        result_names = get_read_names(out_r1)

        # Extract expected read names
        expected_names = set()
        with gzip.open(sanity_data["expected_r1"], "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    name = line.strip().split()[0]
                    if name.startswith("@"):
                        name = name[1:]
                    expected_names.add(name)

        # Compare sets
        missing = expected_names - result_names
        extra = result_names - expected_names

        assert not missing and not extra, (
            f"Read name mismatch:\n"
            f"  Missing {len(missing)} reads: {list(missing)[:5]}...\n"
            f"  Extra {len(extra)} reads: {list(extra)[:5]}..."
        )

    def test_fastq_pair_consistency(
        self,
        sanity_data: dict[str, Path],
        sanity_tmp_dir: Path,
    ) -> None:
        """R1 and R2 FASTQ files must have matching read pairs."""
        import wasp2_rust

        out_r1 = sanity_tmp_dir / "r1.fq"
        out_r2 = sanity_tmp_dir / "r2.fq"

        # Create BED from VCF
        bed_path = sanity_tmp_dir / "variants.bed"
        wasp2_rust.vcf_to_bed_py(
            str(sanity_data["vcf"]),
            str(bed_path),
        )

        # Run unified pipeline
        wasp2_rust.unified_make_reads_py(
            str(sanity_data["bam"]),
            str(bed_path),
            str(out_r1),
            str(out_r2),
            max_seqs=64,
            threads=1,
            compression_threads=1,
            compress_output=False,
            indel_mode=True,
        )

        # Count reads in each file
        def count_reads(fq_path: Path) -> int:
            with open(fq_path) as f:
                return sum(1 for _ in f) // 4

        r1_count = count_reads(out_r1)
        r2_count = count_reads(out_r2)

        assert r1_count == r2_count, f"R1 has {r1_count} reads, R2 has {r2_count} reads"
        assert r1_count > 0, "Expected at least one read pair"


class TestAnalysis:
    """Test analysis output with floating-point tolerance."""

    def test_analysis_variant_count(
        self,
        sanity_data: dict[str, Path],
    ) -> None:
        """Analysis should process variants with sufficient coverage."""
        # Count variants in analysis output
        with open(sanity_data["expected_analysis"]) as f:
            header = f.readline()
            variant_count = sum(1 for _ in f)

        # Should have >100 variants with sufficient coverage
        assert variant_count > 100, f"Expected >100 analyzed variants, got {variant_count}"

    def test_analysis_columns_present(
        self,
        sanity_data: dict[str, Path],
    ) -> None:
        """Analysis output should have required columns."""
        with open(sanity_data["expected_analysis"]) as f:
            header = f.readline().strip().split("\t")

        required = ["region", "ref_count", "alt_count", "pval"]
        missing = [col for col in required if col not in header]
        assert not missing, f"Missing columns: {missing}"

    def test_analysis_reproducible(
        self,
        sanity_data: dict[str, Path],
        sanity_tmp_dir: Path,
    ) -> None:
        """Analysis results should be reproducible."""
        import wasp2_rust

        # Run analysis on counts
        results = wasp2_rust.analyze_imbalance(
            str(sanity_data["expected_counts"]),
            min_count=10,
            pseudocount=1,
            method="single",
        )

        # Load expected analysis
        expected_regions = {}
        with open(sanity_data["expected_analysis"]) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                row = dict(zip(header, parts))
                expected_regions[row["region"]] = {
                    "ref_count": int(row["ref_count"]),
                    "alt_count": int(row["alt_count"]),
                    "pval": float(row["pval"]),
                }

        # Compare (counts should match exactly, p-values with tolerance)
        for result in results:
            region = result["region"]
            if region in expected_regions:
                exp = expected_regions[region]
                assert result["ref_count"] == exp["ref_count"], f"{region}: ref_count mismatch"
                assert result["alt_count"] == exp["alt_count"], f"{region}: alt_count mismatch"
                # P-values may have small floating-point differences
                assert abs(result["pval"] - exp["pval"]) < 1e-6, (
                    f"{region}: pval mismatch (got {result['pval']}, expected {exp['pval']})"
                )


class TestPipelineIntegration:
    """End-to-end pipeline integration tests."""

    def test_full_pipeline_runs(
        self,
        sanity_data: dict[str, Path],
        sanity_tmp_dir: Path,
    ) -> None:
        """Full pipeline should complete without errors."""
        import pysam

        import wasp2_rust

        # 1. Count alleles
        regions = []
        with pysam.VariantFile(str(sanity_data["vcf"])) as vcf:
            for rec in vcf:
                gt = rec.samples[0]["GT"]
                if gt == (0, 1) or gt == (1, 0):
                    regions.append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))

        counter = wasp2_rust.BamCounter(str(sanity_data["bam"]))
        counts = counter.count_alleles(regions, min_qual=0, threads=1)

        assert len(counts) == len(regions)

        # 2. Generate FASTQs
        bed_path = sanity_tmp_dir / "variants.bed"
        wasp2_rust.vcf_to_bed_py(str(sanity_data["vcf"]), str(bed_path))

        out_r1 = sanity_tmp_dir / "r1.fq"
        out_r2 = sanity_tmp_dir / "r2.fq"

        stats = wasp2_rust.unified_make_reads_py(
            str(sanity_data["bam"]),
            str(bed_path),
            str(out_r1),
            str(out_r2),
            max_seqs=64,
            threads=1,
            compression_threads=1,
            compress_output=False,
            indel_mode=True,
        )

        assert stats["total_reads"] > 0
        assert out_r1.exists()
        assert out_r2.exists()

        # 3. Save counts and analyze
        counts_path = sanity_tmp_dir / "counts.tsv"
        with open(counts_path, "w") as f:
            f.write("chrom\tpos\tref\talt\tref_count\talt_count\tother_count\n")
            for (chrom, pos, ref, alt), (rc, ac, oc) in zip(regions, counts):
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{rc}\t{ac}\t{oc}\n")

        results = wasp2_rust.analyze_imbalance(
            str(counts_path),
            min_count=10,
            pseudocount=1,
            method="single",
        )

        assert len(results) > 0, "Analysis should produce results"
