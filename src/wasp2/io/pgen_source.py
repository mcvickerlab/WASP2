"""
PGEN variant source for WASP2.

This module provides a VariantSource implementation for reading PLINK2 PGEN files
using the pgenlib library for efficient genotype access.
"""

import logging
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd

from .variant_source import (
    Genotype,
    Variant,
    VariantGenotype,
    VariantSource,
)

logger = logging.getLogger(__name__)

# Try to import pgenlib - graceful degradation if not available
try:
    import pgenlib

    PGENLIB_AVAILABLE = True
except ImportError:
    PGENLIB_AVAILABLE = False
    logger.debug("pgenlib not available - PGEN functionality will be limited")


@VariantSource.register("pgen")
class PGENSource(VariantSource):
    """PGEN file reader for WASP2.

    Reads PLINK2 PGEN format files using pgenlib for efficient genotype access.
    Automatically locates companion .pvar and .psam files.

    Supports:
    - Multiallelic variants
    - Missing genotypes
    - Heterozygous filtering
    - Region queries
    - BED export

    Args:
        path: Path to .pgen file (or prefix without extension)
        **kwargs: Additional arguments (reserved for future use)

    Raises:
        ImportError: If pgenlib is not installed
        FileNotFoundError: If .pgen, .pvar, or .psam files are missing
        RuntimeError: If PGEN file cannot be opened

    Example:
        >>> source = PGENSource("data/genotypes.pgen")
        >>> for vg in source.iter_variants(het_only=True):
        ...     print(f"{vg.variant.chrom}:{vg.variant.pos}")
    """

    def __init__(self, path: Path, **kwargs):
        """Initialize PGEN source.

        Args:
            path: Path to .pgen file
            **kwargs: Additional arguments (reserved)
        """
        if not PGENLIB_AVAILABLE:
            raise ImportError(
                "pgenlib is required for PGEN support. Install with: pip install pgenlib"
            )

        # Store path and auto-detect companion files
        self.path = Path(path)
        self._detect_companion_files()

        # Read PSAM and PVAR metadata
        self._psam_df = self._read_psam()
        self._pvar_df = self._read_pvar()

        # Initialize pgenlib reader with multiallelic support
        self._reader = self._open_pgen_reader()

    def _detect_companion_files(self):
        """Detect .pvar and .psam files from .pgen path."""
        # If path has .pgen extension, use it directly
        if self.path.suffix == ".pgen":
            pgen_path = self.path
            prefix = self.path.with_suffix("")
        else:
            # Assume path is a prefix
            prefix = self.path
            pgen_path = prefix.with_suffix(".pgen")

        # Set companion file paths
        self.pgen_path = pgen_path
        self.pvar_path = prefix.with_suffix(".pvar")
        self.psam_path = prefix.with_suffix(".psam")

        # Validate all files exist
        if not self.pgen_path.exists():
            raise FileNotFoundError(f"PGEN file not found: {self.pgen_path}")
        if not self.pvar_path.exists():
            raise FileNotFoundError(f"PVAR file not found: {self.pvar_path}")
        if not self.psam_path.exists():
            raise FileNotFoundError(f"PSAM file not found: {self.psam_path}")

    def _read_psam(self) -> pd.DataFrame:
        """Read PSAM file with sample information.

        Returns:
            DataFrame with sample metadata
        """
        # PSAM files may have '#' prefix on header line
        with open(self.psam_path) as f:
            first_line = f.readline().strip()
            has_header = first_line.startswith("#")

        if has_header:
            # Read with header, removing '#' prefix
            df = pd.read_csv(self.psam_path, sep="\t", dtype=str)
            df.columns = [col.lstrip("#") for col in df.columns]
        else:
            # Use default PLINK2 column names
            df = pd.read_csv(self.psam_path, sep="\t", names=["FID", "IID"], dtype=str)

        return df

    def _read_pvar(self) -> pd.DataFrame:
        """Read PVAR file with variant information.

        Returns:
            DataFrame with variant metadata
        """
        # PVAR files have ## comments and optional # header
        # Skip ## lines, but keep # header line
        with open(self.pvar_path) as f:
            lines = f.readlines()

        # Find first non-## line
        data_start = 0
        for i, line in enumerate(lines):
            if not line.startswith("##"):
                data_start = i
                break

        # Check if first data line is header (starts with #CHROM or #)
        has_header = lines[data_start].startswith("#")

        if has_header:
            # Read from data_start, treating first line as header
            df = pd.read_csv(
                self.pvar_path,
                sep="\t",
                skiprows=data_start,
                dtype={"CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str},
            )
            df.columns = [col.lstrip("#") for col in df.columns]
        else:
            # No header - use standard column names
            df = pd.read_csv(
                self.pvar_path,
                sep="\t",
                skiprows=data_start,
                names=["CHROM", "POS", "ID", "REF", "ALT"],
                dtype={"CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str},
            )

        # Normalize chromosome names to include 'chr' prefix for consistency
        # plink2 strips 'chr' prefix by default, but we want consistent output
        df["CHROM"] = df["CHROM"].apply(self._normalize_chrom_name)

        return df

    def _normalize_chrom_name(self, chrom: str) -> str:
        """Normalize chromosome name to include 'chr' prefix.

        Args:
            chrom: Chromosome name (e.g., '1', 'chr1', 'X')

        Returns:
            Normalized chromosome name with 'chr' prefix
        """
        chrom = str(chrom)
        # Already has chr prefix
        if chrom.lower().startswith("chr"):
            return chrom
        # Add chr prefix for numeric chromosomes
        if chrom.isdigit() or chrom in ("X", "Y", "M", "MT"):
            return f"chr{chrom}"
        return chrom

    def _open_pgen_reader(self):
        """Open pgenlib reader with multiallelic support.

        Returns:
            pgenlib.PgenReader instance
        """
        # Calculate allele counts for multiallelic support
        # Count commas in ALT field + 2 (REF + ALT alleles)
        allele_counts = self._pvar_df["ALT"].str.count(",") + 2

        # Create allele index offsets for pgenlib
        allele_idx_offsets = np.zeros(len(self._pvar_df) + 1, dtype=np.uintp)
        allele_idx_offsets[1:] = np.cumsum(allele_counts)

        try:
            # pgenlib expects bytes for filename
            reader = pgenlib.PgenReader(
                bytes(str(self.pgen_path), "utf-8"), allele_idx_offsets=allele_idx_offsets
            )
            return reader
        except Exception as e:
            raise RuntimeError(f"Failed to open PGEN file: {e}") from e

    @property
    def samples(self) -> list[str]:
        """Get list of sample IDs.

        Returns:
            List of sample IDs from PSAM file
        """
        # Try common sample ID columns
        for col in ["IID", "ID", "SAMPLE"]:
            if col in self._psam_df.columns:
                return list(self._psam_df[col].tolist())

        # Fallback to first column
        return list(self._psam_df.iloc[:, 0].tolist())

    @property
    def variant_count(self) -> int:
        """Get total number of variants.

        Returns:
            Number of variants in PGEN file
        """
        return int(self._reader.get_variant_ct())

    @property
    def sample_count(self) -> int:
        """Get total number of samples.

        Returns:
            Number of samples in PGEN file
        """
        return int(self._reader.get_raw_sample_ct())

    def iter_variants(
        self, samples: list[str] | None = None, het_only: bool = False
    ) -> Iterator[VariantGenotype]:
        """Iterate over variants with optional filtering.

        Args:
            samples: Optional list of sample IDs to include. If None, use first sample.
            het_only: If True, only yield heterozygous variants

        Yields:
            VariantGenotype objects for each variant/sample combination
        """
        # Determine which samples to process
        if samples is None:
            # Default to first sample
            sample_indices = [0]
            sample_ids = [self.samples[0]]
        else:
            sample_indices = [self.get_sample_idx(s) for s in samples]
            sample_ids = samples

        # Iterate through all variants
        for variant_idx in range(self.variant_count):
            variant_row = self._pvar_df.iloc[variant_idx]

            # Create Variant object
            variant = Variant(
                chrom=str(variant_row["CHROM"]),
                pos=int(variant_row["POS"]),
                ref=str(variant_row["REF"]),
                alt=str(variant_row["ALT"]),
                id=str(variant_row["ID"]) if "ID" in variant_row else None,
            )

            # Read genotypes for each requested sample
            for sample_idx, _sample_id in zip(sample_indices, sample_ids):
                # Set sample subset for this sample
                sample_subset = np.array([sample_idx], dtype=np.uint32)
                self._reader.change_sample_subset(sample_subset)

                # Read alleles for this variant
                allele_buf = np.zeros(2, dtype=np.int32)
                self._reader.read_alleles(variant_idx, allele_buf)

                # Parse genotype
                genotype, allele1, allele2 = self._parse_alleles(allele_buf, variant_row)

                # Apply het_only filter
                if het_only and genotype != Genotype.HET:
                    continue

                # Yield VariantGenotype
                yield VariantGenotype(
                    variant=variant, genotype=genotype, allele1=allele1, allele2=allele2
                )

    def get_genotype(self, sample: str, chrom: str, pos: int) -> Genotype:
        """Get genotype for a specific sample at a genomic position.

        Args:
            sample: Sample ID
            chrom: Chromosome name
            pos: 1-based genomic position

        Returns:
            Genotype enum value

        Raises:
            ValueError: If sample not found or position has no variant
        """
        # Find sample index
        sample_idx = self.get_sample_idx(sample)

        # Normalize chromosome for comparison (handle both str and int)
        chrom_normalized = self._normalize_chrom(chrom)

        # Find variant by chrom/pos
        mask = (self._pvar_df["CHROM"] == chrom_normalized) & (self._pvar_df["POS"] == pos)
        matching_variants = self._pvar_df[mask]

        if len(matching_variants) == 0:
            raise ValueError(f"No variant found at {chrom}:{pos}")

        variant_idx = matching_variants.index[0]
        variant_row = matching_variants.iloc[0]

        # Set sample subset and read genotype
        sample_subset = np.array([sample_idx], dtype=np.uint32)
        self._reader.change_sample_subset(sample_subset)

        allele_buf = np.zeros(2, dtype=np.int32)
        self._reader.read_alleles(variant_idx, allele_buf)

        # Parse and return genotype
        genotype, _, _ = self._parse_alleles(allele_buf, variant_row)
        return genotype

    def query_region(
        self, chrom: str, start: int, end: int, samples: list[str] | None = None
    ) -> Iterator[VariantGenotype]:
        """Query variants in a genomic region.

        Uses 1-based inclusive coordinates.

        Args:
            chrom: Chromosome name
            start: 1-based start position (inclusive)
            end: 1-based end position (inclusive)
            samples: Optional list of sample IDs to include

        Yields:
            VariantGenotype objects in the region
        """
        # Normalize chromosome for comparison (handle both str and int)
        chrom_normalized = self._normalize_chrom(chrom)

        # Filter PVAR by region
        mask = (
            (self._pvar_df["CHROM"] == chrom_normalized)
            & (self._pvar_df["POS"] >= start)
            & (self._pvar_df["POS"] <= end)
        )
        region_variants = self._pvar_df[mask]

        # Determine samples
        if samples is None:
            sample_indices = [0]
            sample_ids = [self.samples[0]]
        else:
            sample_indices = [self.get_sample_idx(s) for s in samples]
            sample_ids = samples

        # Iterate through variants in region
        for idx in region_variants.index:
            variant_row = self._pvar_df.loc[idx]

            variant = Variant(
                chrom=str(variant_row["CHROM"]),
                pos=int(variant_row["POS"]),
                ref=str(variant_row["REF"]),
                alt=str(variant_row["ALT"]),
                id=str(variant_row["ID"]) if "ID" in variant_row else None,
            )

            # Read genotypes for requested samples
            for sample_idx, _sample_id in zip(sample_indices, sample_ids):
                sample_subset = np.array([sample_idx], dtype=np.uint32)
                self._reader.change_sample_subset(sample_subset)

                allele_buf = np.zeros(2, dtype=np.int32)
                self._reader.read_alleles(idx, allele_buf)

                genotype, allele1, allele2 = self._parse_alleles(allele_buf, variant_row)

                yield VariantGenotype(
                    variant=variant, genotype=genotype, allele1=allele1, allele2=allele2
                )

    def to_bed(
        self,
        output: Path,
        samples: list[str] | None = None,
        het_only: bool = True,
        include_genotypes: bool = True,
    ) -> Path:
        """Export variants to BED format file.

        BED format uses 0-based start, 1-based end coordinates.

        Args:
            output: Output BED file path
            samples: Optional list of sample IDs to include
            het_only: If True, only export heterozygous variants
            include_genotypes: If True, include genotype column

        Returns:
            Path to the created BED file
        """
        output_path = Path(output)

        with open(output_path, "w") as f:
            for vg in self.iter_variants(samples=samples, het_only=het_only):
                # Write BED line: chrom, start (0-based), end (1-based), ref, alt
                line = vg.variant.to_bed_line()

                # Add genotype if requested
                if include_genotypes:
                    gt_str = self._genotype_to_string(vg.genotype)
                    line += f"\t{gt_str}"

                f.write(line + "\n")

        return output_path

    def _normalize_chrom(self, chrom: str) -> str:
        """Normalize chromosome value for queries.

        Since we normalize PVAR chromosomes to have 'chr' prefix,
        we need to normalize query chromosomes the same way.

        Args:
            chrom: Chromosome name (str or int-like)

        Returns:
            Normalized chromosome value with 'chr' prefix
        """
        return self._normalize_chrom_name(str(chrom))

    def _parse_alleles(
        self, allele_buf: np.ndarray, variant_row
    ) -> tuple[Genotype, str | None, str | None]:
        """Convert allele buffer to Genotype and allele sequences.

        Args:
            allele_buf: Array with two allele indices
            variant_row: PVAR row for this variant

        Returns:
            Tuple of (Genotype, allele1_seq, allele2_seq)
        """
        allele1_idx = allele_buf[0]
        allele2_idx = allele_buf[1]

        # Check for missing genotype (-9 in pgenlib)
        if allele1_idx < 0 or allele2_idx < 0:
            return Genotype.MISSING, None, None

        # Get allele sequences
        allele1_seq = self._allele_idx_to_base(allele1_idx, variant_row)
        allele2_seq = self._allele_idx_to_base(allele2_idx, variant_row)

        # Classify genotype
        if allele1_idx == allele2_idx:
            if allele1_idx == 0:
                return Genotype.HOM_REF, allele1_seq, allele2_seq
            else:
                return Genotype.HOM_ALT, allele1_seq, allele2_seq
        else:
            return Genotype.HET, allele1_seq, allele2_seq

    def _allele_idx_to_base(self, idx: int, variant_row) -> str:
        """Convert allele index to base sequence.

        Args:
            idx: Allele index (0=REF, 1+=ALT)
            variant_row: PVAR row for this variant

        Returns:
            Allele sequence string
        """
        if idx == 0:
            return str(variant_row["REF"])
        else:
            # ALT may be comma-separated for multiallelic
            alt_alleles = str(variant_row["ALT"]).split(",")
            alt_idx = idx - 1
            if alt_idx < len(alt_alleles):
                return alt_alleles[alt_idx]
            else:
                # Should not happen with correct allele_idx_offsets
                logger.warning(f"Invalid ALT index {alt_idx} for variant")
                return "."

    def _genotype_to_string(self, genotype: Genotype) -> str:
        """Convert Genotype enum to string representation.

        Args:
            genotype: Genotype enum value

        Returns:
            String representation (e.g., "0/1", "1/1")
        """
        if genotype == Genotype.HOM_REF:
            return "0/0"
        elif genotype == Genotype.HET:
            return "0/1"
        elif genotype == Genotype.HOM_ALT:
            return "1/1"
        else:
            return "./."

    def close(self):
        """Close the PGEN reader and release resources."""
        if hasattr(self, "_reader") and self._reader is not None:
            self._reader.close()
            self._reader = None
