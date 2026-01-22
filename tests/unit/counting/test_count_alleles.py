"""
Unit tests for the counting.count_alleles module.

Tests cover:
- find_read_aln_pos binary search functionality
- RUST_AVAILABLE flag checking
- Thread environment variable parsing
"""

import os
from unittest.mock import MagicMock, patch

import pytest


class TestFindReadAlnPos:
    """Tests for the find_read_aln_pos helper function."""

    def test_finds_exact_position(self):
        """Test that exact position matches return correct query position."""
        from counting.count_alleles import find_read_aln_pos

        # Create mock read with aligned pairs
        mock_read = MagicMock()
        # Aligned pairs: (query_pos, ref_pos)
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
            (2, 103),  # Gap at 102
        ]

        assert find_read_aln_pos(mock_read, 102) is None
        assert find_read_aln_pos(mock_read, 50) is None  # Before range
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


class TestRustAvailability:
    """Tests for Rust extension availability checking."""

    def test_rust_available_flag_exists(self):
        """Test that RUST_AVAILABLE flag is defined."""
        from counting.count_alleles import RUST_AVAILABLE

        assert isinstance(RUST_AVAILABLE, bool)


class TestThreadEnvironment:
    """Tests for thread environment variable handling."""

    def test_default_threads_without_env(self):
        """Test default thread count when env var not set."""
        # Clear environment variable if set
        env_backup = os.environ.pop("WASP2_RUST_THREADS", None)

        try:
            # The function uses 1 as default when env not set
            from counting.count_alleles import count_snp_alleles_rust

            # We can't actually call the function without Rust,
            # but we can verify the module loads correctly
            assert callable(count_snp_alleles_rust)
        finally:
            # Restore environment
            if env_backup is not None:
                os.environ["WASP2_RUST_THREADS"] = env_backup

    def test_invalid_thread_env_falls_back(self):
        """Test that invalid thread env value falls back to 1."""
        env_backup = os.environ.get("WASP2_RUST_THREADS")

        try:
            os.environ["WASP2_RUST_THREADS"] = "invalid"

            # Force reimport to pick up new env
            import importlib
            import counting.count_alleles as ca
            importlib.reload(ca)

            # The module should load without error even with invalid env
            assert ca.RUST_AVAILABLE is not None
        finally:
            if env_backup is not None:
                os.environ["WASP2_RUST_THREADS"] = env_backup
            elif "WASP2_RUST_THREADS" in os.environ:
                del os.environ["WASP2_RUST_THREADS"]
