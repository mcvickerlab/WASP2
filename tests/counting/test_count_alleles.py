import polars as pl
import pytest

from counting import count_alleles


def _variants() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": ["chr1"],
            "pos": [10],
            "ref": ["A"],
            "alt": ["G"],
        }
    ).with_columns(pl.col("chrom").cast(pl.Categorical))


@pytest.mark.unit
def test_bulk_counting_fails_closed_on_chromosome_error(monkeypatch):
    monkeypatch.setattr(count_alleles, "RUST_AVAILABLE", True)

    def fail(*args, **kwargs):
        raise RuntimeError("counter failed")

    monkeypatch.setattr(count_alleles, "count_snp_alleles_rust", fail)
    with pytest.raises(RuntimeError, match="counter failed"):
        count_alleles.make_count_df("sample.bam", _variants())


@pytest.mark.unit
def test_bulk_counting_preserves_counts_above_uint16(monkeypatch):
    monkeypatch.setattr(count_alleles, "RUST_AVAILABLE", True)
    monkeypatch.setattr(
        count_alleles,
        "count_snp_alleles_rust",
        lambda *args, **kwargs: [("chr1", 10, 70_000, 80_000, 90_000)],
    )
    result = count_alleles.make_count_df("sample.bam", _variants())
    assert result.select("ref_count", "alt_count", "other_count").row(0) == (
        70_000,
        80_000,
        90_000,
    )
    assert result.schema["ref_count"] == pl.UInt32
    assert result.schema["alt_count"] == pl.UInt32
    assert result.schema["other_count"] == pl.UInt32
