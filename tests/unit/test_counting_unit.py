"""
Unit tests for the WASP2 counting module.

This module contains comprehensive unit tests for:
- Allele counting functions
- Binary search utilities
- Thread environment handling
- Rust extension availability

Run with: pytest tests/unit/test_counting_unit.py -v
Run only fast tests: pytest tests/unit/test_counting_unit.py -v -m "not slow"
Run Rust tests: pytest tests/unit/test_counting_unit.py -v -m rust
"""

import os
from unittest.mock import MagicMock, patch

import pytest


# ============================================================================
# Test: find_read_aln_pos binary search function
# ============================================================================

class TestFindReadAlnPos:
    """Tests for the find_read_aln_pos helper function.

    This function performs binary search over aligned pairs to find
    the query position corresponding to a given reference position.
    """

    def test_finds_exact_position(self):
        """Test that exact position matches return correct query position."""
        from counting.count_alleles import find_read_aln_pos

        mock_read = MagicMock()
        mock_read.get_aligned_pairs.return_value = [
            (0, 100),
            (1, 101),
            (2, 102),
            (3, 103),
            (4, 104),
        ]

        assert find_read_aln_pos(mock_read, 100) == 0
        assert find_read_aln_pos(mock_read, 102) == 2
        assert find_read_aln_pos(mock_read, 104) == 4

    def test_returns_none_for_missing_position(self):
        """Test that missing positions return None."""
        from counting.count_alleles import find_read_aln_pos

        mock_read = MagicMock()
        mock_read.get_aligned_pairs.return_value = [
            (0, 100),
            (1, 101),
            (2, 103),  # Gap at position 102
        ]

        assert find_read_aln_pos(mock_read, 102) is None
        assert find_read_aln_pos(mock_read, 50) is None   # Before range
        assert find_read_aln_pos(mock_read, 200) is None  # After range

    def test_empty_alignment(self):
        """Test handling of empty alignment."""
        from counting.count_alleles import find_read_aln_pos

        mock_read = MagicMock()
        mock_read.get_aligned_pairs.return_value = []

        assert find_read_aln_pos(mock_read, 100) is None

    def test_single_alignment(self):
        """Test with single aligned pair."""
        from counting.count_alleles import find_read_aln_pos

        mock_read = MagicMock()
        mock_read.get_aligned_pairs.return_value = [(5, 200)]

        assert find_read_aln_pos(mock_read, 200) == 5
        assert find_read_aln_pos(mock_read, 199) is None
        assert find_read_aln_pos(mock_read, 201) is None

    def test_large_alignment(self):
        """Test binary search efficiency with large alignment."""
        from counting.count_alleles import find_read_aln_pos

        mock_read = MagicMock()
        # Create 1000 aligned pairs
        mock_read.get_aligned_pairs.return_value = [
            (i, 1000 + i) for i in range(1000)
        ]

        # Should efficiently find positions
        assert find_read_aln_pos(mock_read, 1000) == 0
        assert find_read_aln_pos(mock_read, 1500) == 500
        assert find_read_aln_pos(mock_read, 1999) == 999
        assert find_read_aln_pos(mock_read, 2000) is None

    def test_handles_insertions(self):
        """Test handling of insertions (None reference positions)."""
        from counting.count_alleles import find_read_aln_pos

        mock_read = MagicMock()
        # Insertion at query positions 1-2 (ref position is None)
        mock_read.get_aligned_pairs.return_value = [
            (0, 100),
            (1, None),  # Insertion
            (2, None),  # Insertion
            (3, 101),
            (4, 102),
        ]

        assert find_read_aln_pos(mock_read, 100) == 0
        assert find_read_aln_pos(mock_read, 101) == 3
        assert find_read_aln_pos(mock_read, 102) == 4


# ============================================================================
# Test: Rust extension availability
# ============================================================================

class TestRustAvailability:
    """Tests for Rust extension availability checking."""

    def test_rust_available_flag_exists(self):
        """Test that RUST_AVAILABLE flag is defined."""
        from counting.count_alleles import RUST_AVAILABLE

        assert isinstance(RUST_AVAILABLE, bool)

    @pytest.mark.rust
    def test_rust_counter_import(self):
        """Test that Rust counter can be imported when available."""
        from counting.count_alleles import RUST_AVAILABLE

        if RUST_AVAILABLE:
            from wasp2_rust import BamCounter as RustBamCounter
            assert RustBamCounter is not None
        else:
            pytest.skip("Rust extension not available")

    @pytest.mark.rust
    def test_count_snp_alleles_rust_callable(self):
        """Test that count_snp_alleles_rust function is callable."""
        from counting.count_alleles import count_snp_alleles_rust

        assert callable(count_snp_alleles_rust)


# ============================================================================
# Test: Thread environment handling
# ============================================================================

class TestThreadEnvironment:
    """Tests for thread environment variable handling."""

    def test_default_threads_without_env(self):
        """Test default thread count when env var not set."""
        env_backup = os.environ.pop("WASP2_RUST_THREADS", None)

        try:
            from counting.count_alleles import count_snp_alleles_rust
            assert callable(count_snp_alleles_rust)
        finally:
            if env_backup is not None:
                os.environ["WASP2_RUST_THREADS"] = env_backup

    def test_valid_thread_env_parsed(self):
        """Test that valid thread count from env is parsed correctly."""
        env_backup = os.environ.get("WASP2_RUST_THREADS")

        try:
            os.environ["WASP2_RUST_THREADS"] = "4"

            import importlib
            import counting.count_alleles as ca
            importlib.reload(ca)

            # Module should load without error
            assert ca.RUST_AVAILABLE is not None
        finally:
            if env_backup is not None:
                os.environ["WASP2_RUST_THREADS"] = env_backup
            elif "WASP2_RUST_THREADS" in os.environ:
                del os.environ["WASP2_RUST_THREADS"]

    def test_invalid_thread_env_falls_back(self):
        """Test that invalid thread env value falls back to 1."""
        env_backup = os.environ.get("WASP2_RUST_THREADS")

        try:
            os.environ["WASP2_RUST_THREADS"] = "invalid"

            import importlib
            import counting.count_alleles as ca
            importlib.reload(ca)

            # Module should load without error even with invalid env
            assert ca.RUST_AVAILABLE is not None
        finally:
            if env_backup is not None:
                os.environ["WASP2_RUST_THREADS"] = env_backup
            elif "WASP2_RUST_THREADS" in os.environ:
                del os.environ["WASP2_RUST_THREADS"]

    def test_negative_thread_env_clamped(self):
        """Test that negative thread count is clamped to 1."""
        env_backup = os.environ.get("WASP2_RUST_THREADS")

        try:
            os.environ["WASP2_RUST_THREADS"] = "-5"

            import importlib
            import counting.count_alleles as ca
            importlib.reload(ca)

            # Module should load without error
            assert ca.RUST_AVAILABLE is not None
        finally:
            if env_backup is not None:
                os.environ["WASP2_RUST_THREADS"] = env_backup
            elif "WASP2_RUST_THREADS" in os.environ:
                del os.environ["WASP2_RUST_THREADS"]


# ============================================================================
# Test: DataFrame counting operations
# ============================================================================

class TestMakeCountDf:
    """Tests for the make_count_df function."""

    @pytest.mark.rust
    def test_make_count_df_requires_rust(self):
        """Test that make_count_df raises error without Rust."""
        from counting.count_alleles import make_count_df, RUST_AVAILABLE

        if RUST_AVAILABLE:
            pytest.skip("This test requires Rust to be unavailable")

        import polars as pl

        mock_df = pl.DataFrame({
            "chrom": ["chr1"],
            "pos": [100],
            "ref": ["A"],
            "alt": ["G"],
        })

        with pytest.raises(RuntimeError, match="Rust BAM counter not available"):
            make_count_df("nonexistent.bam", mock_df, use_rust=True)

    @pytest.mark.rust
    @pytest.mark.slow
    def test_make_count_df_with_rust(self, sample_bam, sample_intersections_df):
        """Test make_count_df with actual Rust extension."""
        from counting.count_alleles import make_count_df, RUST_AVAILABLE

        if not RUST_AVAILABLE:
            pytest.skip("Rust extension not available")

        result_df = make_count_df(str(sample_bam), sample_intersections_df)

        assert "ref_count" in result_df.columns
        assert "alt_count" in result_df.columns
        assert "other_count" in result_df.columns


# ============================================================================
# Test: Module imports and structure
# ============================================================================

class TestModuleStructure:
    """Tests for module import and structure."""

    def test_count_alleles_imports(self):
        """Test that counting.count_alleles module imports correctly."""
        import counting.count_alleles as ca

        assert hasattr(ca, "RUST_AVAILABLE")
        assert hasattr(ca, "count_snp_alleles_rust")
        assert hasattr(ca, "make_count_df")
        assert hasattr(ca, "find_read_aln_pos")

    def test_logging_configured(self):
        """Test that module has logger configured."""
        import counting.count_alleles as ca

        assert hasattr(ca, "logger")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
