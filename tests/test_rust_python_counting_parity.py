"""Bulk SNP counting parity against archived Python commit
41354fc087404ba1fa99542803ee9cb9b7e0e6bf.

The fixture is the verbatim archived src/counting/count_alleles.py. Only its
count_snp_alleles function is loaded, avoiding unrelated legacy imports.
"""

import ast
from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[1]
SHARED_DATA = ROOT / "tests" / "shared_data"
SANITY_DATA = ROOT / "tests" / "sanity" / "data"


# ---------------------------------------------------------------------------
# Unmodified archived Python counting reference
# ---------------------------------------------------------------------------
def python_count_snp_alleles(bam_path: str, chrom: str, snp_list: list[tuple]):
    """Run the unmodified archived function, not a translation of Rust."""
    source = (ROOT / "tests/counting/fixtures/archive_count_alleles.py.txt").read_text()
    tree = ast.parse(source)
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "count_snp_alleles"
    )
    namespace = {}
    exec(
        compile(ast.Module(body=[function], type_ignores=[]), "archived_python", "exec"), namespace
    )
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return namespace["count_snp_alleles"](bam, chrom, snp_list)


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

        Compare the available HG00731 data to the archived algorithm.
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


def _write_counting_bam(tmp_path, scenario):
    bam_path = tmp_path / "archive_priority.bam"
    # Records are coordinate sorted. The same name on chr2 must remain eligible.
    records = [("chr1", 100, "shared", "A" * 20, [(0, 20)], 65)]
    if scenario == "mates":
        records.append(("chr1", 200, "shared", "G" * 20, [(0, 20)], 129))
    elif scenario == "spanning":
        records = [("chr1", 100, "shared", "A" * 120, [(0, 120)], 0)]
    elif scenario in {"deletion", "skip"}:
        gap = 2 if scenario == "deletion" else 3
        records = [("chr1", 99, "shared", "A" * 10, [(0, 1), (gap, 100), (0, 9)], 0)]
    elif scenario == "repeated":
        records.append(("chr1", 100, "shared", "G" * 20, [(0, 20)], 0))
        records.append(("chr1", 200, "shared", "G" * 20, [(0, 20)], 0))
    elif scenario == "other":
        records = [("chr1", 100, "shared", "N" + "A" * 19, [(0, 20)], 0)]
    records.append(("chr2", 100, "shared", "G" * 20, [(0, 20)], 0))
    with pysam.AlignmentFile(
        str(bam_path),
        "wb",
        header={
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}],
        },
    ) as bam:
        for chrom, start, name, sequence, cigar, flag in records:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = sequence
            read.flag = flag
            read.reference_id = 0 if chrom == "chr1" else 1
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = cigar
            read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            if flag & 1:
                read.next_reference_id = 0
                read.next_reference_start = 200 if flag == 65 else 100
            bam.write(read)
    pysam.index(str(bam_path))
    return str(bam_path)


@pytest.mark.parametrize("scenario", ["mates", "spanning", "deletion", "skip", "repeated", "other"])
@pytest.mark.parametrize("positions", [[101, 201], [201, 101]])
@pytest.mark.parametrize("threads", [1, 2])
def test_archived_snp_priority(tmp_path, scenario, positions, threads):
    import wasp2_rust

    bam_path = _write_counting_bam(tmp_path, scenario)
    snps = [(pos, "A", "G") for pos in positions]
    expected = [row[2:] for row in python_count_snp_alleles(bam_path, "chr1", snps)]
    expected += [row[2:] for row in python_count_snp_alleles(bam_path, "chr2", [(101, "A", "G")])]
    regions = [("chr1", *snp) for snp in snps] + [("chr2", 101, "A", "G")]
    observed = wasp2_rust.BamCounter(bam_path).count_alleles(regions, min_qual=0, threads=threads)
    assert observed == expected
    assert observed[-1] == (0, 1, 0), "Read-name deduplication must reset per chromosome"
    if scenario in {"deletion", "skip"} and positions == [101, 201]:
        assert observed[:2] == [(0, 0, 0), (0, 0, 0)]


def test_archived_counts_are_reused_for_overlapping_regions(tmp_path):
    import polars as pl

    from counting.count_alleles import make_count_df

    bam_path = _write_counting_bam(tmp_path, "spanning")
    variants = pl.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "pos": [201, 101, 201, 101],
            "ref": ["A"] * 4,
            "alt": ["G"] * 4,
            "region": ["peak_a", "peak_a", "peak_b", "peak_b"],
        }
    ).with_columns(pl.col("chrom").cast(pl.Categorical))
    expected = {
        row[1]: row[2:]
        for row in python_count_snp_alleles(bam_path, "chr1", [(201, "A", "G"), (101, "A", "G")])
    }
    observed = make_count_df(bam_path, variants)
    assert observed.height == 4
    for row in observed.iter_rows(named=True):
        assert (row["ref_count"], row["alt_count"], row["other_count"]) == expected[row["pos"]]
