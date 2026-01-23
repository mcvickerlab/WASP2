"""
Variant source module for WASP2.

This module provides core data structures and an abstract base class for reading
variant data from different file formats (VCF, PGEN).
"""

from abc import ABC, abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class Genotype(Enum):
    """Genotype encoding for variants.

    Standard VCF-style encoding:
    - HOM_REF: Homozygous reference (0/0)
    - HET: Heterozygous (0/1 or 1/0)
    - HOM_ALT: Homozygous alternate (1/1)
    - MISSING: Missing genotype (./.)
    """

    HOM_REF = 0
    HET = 1
    HOM_ALT = 2
    MISSING = -1


@dataclass(frozen=True, slots=True)
class Variant:
    """Immutable variant data structure.

    Represents a single genomic variant with chromosome, position, and alleles.
    Uses 1-based genomic coordinates (VCF standard).

    Attributes:
        chrom: Chromosome name (e.g., "chr1", "1")
        pos: 1-based genomic position
        ref: Reference allele sequence
        alt: Alternate allele sequence
        id: Optional variant ID (e.g., rsID)
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    id: str | None = None

    @property
    def pos0(self) -> int:
        """Return 0-based position for BED format compatibility.

        Returns:
            0-based position (pos - 1)
        """
        return self.pos - 1

    def to_bed_line(self) -> str:
        """Convert variant to BED format line.

        BED format uses 0-based start, 1-based end coordinates.
        Format: chrom\\tstart\\tend\\tref\\talt

        Returns:
            Tab-separated BED format string
        """
        return f"{self.chrom}\t{self.pos0}\t{self.pos}\t{self.ref}\t{self.alt}"


@dataclass
class VariantGenotype:
    """Variant with genotype information for a specific sample.

    Combines a Variant with genotype data, representing the state
    of this variant in a particular sample.

    Attributes:
        variant: The Variant object
        genotype: Genotype classification (HOM_REF, HET, HOM_ALT, MISSING)
        allele1: Optional first allele sequence
        allele2: Optional second allele sequence
    """

    variant: Variant
    genotype: Genotype
    allele1: str | None = None
    allele2: str | None = None

    @property
    def is_het(self) -> bool:
        """Check if this is a heterozygous genotype.

        Returns:
            True if genotype is HET, False otherwise
        """
        return self.genotype == Genotype.HET


class VariantSource(ABC):
    """Abstract base class for variant file readers with factory pattern.

    VariantSource provides a unified interface for reading variant data from
    different file formats (VCF, PGEN, etc.). It implements a factory pattern
    with automatic format detection and a registry system for format handlers.

    The class supports:
    - Automatic format detection from file extensions
    - Compressed file handling (.gz, .bgz, .zst)
    - Context manager protocol for resource management
    - Iteration over variants with optional filtering
    - Region queries for indexed formats
    - BED format export

    Subclasses must implement:
    - Abstract properties: samples, variant_count, sample_count
    - Abstract methods: iter_variants, get_genotype, query_region, to_bed
    - Optional: close() for cleanup

    Usage:
        # Factory pattern with automatic format detection
        with VariantSource.open("variants.vcf.gz") as source:
            for vg in source.iter_variants(het_only=True):
                print(f"{vg.variant.chrom}:{vg.variant.pos}")

        # Direct subclass instantiation
        from wasp2.io.vcf_source import VCFSource
        source = VCFSource("variants.vcf.gz")
        samples = source.samples
        source.close()

    Registering a new format handler:
        @VariantSource.register("vcf", "bcf")
        class VCFSource(VariantSource):
            def __init__(self, path: str):
                self.path = path
            # ... implement abstract methods
    """

    _registry: dict[str, type] = {}

    @classmethod
    def register(cls, *extensions: str):
        """Decorator to register format handlers for specific file extensions.

        This decorator allows subclasses to register themselves as handlers
        for one or more file extensions. When VariantSource.open() is called,
        the factory will automatically select the appropriate handler based
        on the file extension.

        Args:
            *extensions: Variable number of file extensions (with or without leading dot).
                        Extensions are normalized to lowercase without leading dots.

        Returns:
            Decorator function that registers the subclass and returns it unchanged.

        Example:
            @VariantSource.register("vcf", "bcf")
            class VCFSource(VariantSource):
                pass

            @VariantSource.register(".pgen")
            class PGENSource(VariantSource):
                pass
        """

        def decorator(subclass):
            for ext in extensions:
                cls._registry[ext.lower().lstrip(".")] = subclass
            return subclass

        return decorator

    @classmethod
    def _detect_format(cls, path: Path) -> str:
        """Detect file format from path extension.

        Handles both plain and compressed files. For compressed files
        (.gz, .bgz, .zst), looks at the second-to-last suffix to determine
        the actual format.

        Args:
            path: Path to the variant file

        Returns:
            Format extension as a lowercase string (e.g., "vcf", "pgen")

        Examples:
            >>> VariantSource._detect_format(Path("data.vcf"))
            'vcf'
            >>> VariantSource._detect_format(Path("data.vcf.gz"))
            'vcf'
            >>> VariantSource._detect_format(Path("data.pgen"))
            'pgen'
        """
        suffixes = path.suffixes
        # Compression extensions to skip
        compression_exts = {".gz", ".bgz", ".zst"}

        if not suffixes:
            raise ValueError(f"Cannot detect format: no extension in {path}")

        # If last suffix is compression, use second-to-last
        if len(suffixes) >= 2 and suffixes[-1] in compression_exts:
            return suffixes[-2].lstrip(".").lower()
        else:
            return suffixes[-1].lstrip(".").lower()

    @classmethod
    def open(cls, path: str, **kwargs) -> "VariantSource":
        """Factory method to open a variant file with automatic format detection.

        Automatically detects the file format from the extension and instantiates
        the appropriate handler subclass. Raises descriptive errors if the file
        doesn't exist or the format is not supported.

        Args:
            path: Path to the variant file (str or Path-like)
            **kwargs: Additional arguments passed to the format handler constructor

        Returns:
            Instance of the appropriate VariantSource subclass

        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file format is not supported (no registered handler)

        Examples:
            >>> source = VariantSource.open("data.vcf.gz")
            >>> type(source).__name__
            'VCFSource'

            >>> source = VariantSource.open("data.pgen")
            >>> type(source).__name__
            'PGENSource'
        """
        file_path = Path(path)

        # Check if file exists
        if not file_path.exists():
            raise FileNotFoundError(f"Variant file not found: {path}")

        # Detect format
        format_ext = cls._detect_format(file_path)

        # Look up handler in registry
        if format_ext not in cls._registry:
            supported = ", ".join(sorted(cls._registry.keys()))
            raise ValueError(
                f"Unsupported variant file format: '{format_ext}'. Supported formats: {supported}"
            )

        # Instantiate the appropriate handler
        handler_class = cls._registry[format_ext]
        return handler_class(path, **kwargs)

    @property
    @abstractmethod
    def samples(self) -> list[str]:
        """Get list of sample IDs in the variant file.

        Returns:
            List of sample ID strings in file order
        """
        pass

    @property
    @abstractmethod
    def variant_count(self) -> int:
        """Get total number of variants in the file.

        For some formats, this may require a full file scan if not
        indexed or if the count is not stored in metadata.

        Returns:
            Total number of variants
        """
        pass

    @property
    @abstractmethod
    def sample_count(self) -> int:
        """Get total number of samples in the file.

        Returns:
            Total number of samples
        """
        pass

    @abstractmethod
    def iter_variants(
        self, samples: list[str] | None = None, het_only: bool = False
    ) -> Iterator[VariantGenotype]:
        """Iterate over variants with optional filtering.

        Args:
            samples: Optional list of sample IDs to include. If None, use all samples.
                    For multi-sample iteration, yields one VariantGenotype per sample.
            het_only: If True, only yield heterozygous variants

        Yields:
            VariantGenotype objects for each variant/sample combination

        Example:
            >>> for vg in source.iter_variants(samples=["sample1"], het_only=True):
            ...     print(vg.variant.pos, vg.genotype)
        """
        pass

    @abstractmethod
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
        pass

    @abstractmethod
    def query_region(
        self, chrom: str, start: int, end: int, samples: list[str] | None = None
    ) -> Iterator[VariantGenotype]:
        """Query variants in a genomic region.

        Requires the variant file to be indexed (e.g., .tbi, .csi for VCF).
        Uses 1-based inclusive coordinates.

        Args:
            chrom: Chromosome name
            start: 1-based start position (inclusive)
            end: 1-based end position (inclusive)
            samples: Optional list of sample IDs to include

        Yields:
            VariantGenotype objects in the region

        Raises:
            ValueError: If the file is not indexed or region is invalid
        """
        pass

    @abstractmethod
    def to_bed(
        self,
        output: Path,
        samples: list[str] | None = None,
        het_only: bool = True,
        include_genotypes: bool = True,
    ) -> Path:
        """Export variants to BED format file.

        BED format uses 0-based start, 1-based end coordinates.
        Format depends on include_genotypes:
        - If True: chrom\\tstart\\tend\\tref\\talt\\tgenotype
        - If False: chrom\\tstart\\tend\\tref\\talt

        Args:
            output: Output BED file path
            samples: Optional list of sample IDs to include
            het_only: If True, only export heterozygous variants
            include_genotypes: If True, include genotype column

        Returns:
            Path to the created BED file

        Raises:
            IOError: If file cannot be written
        """
        pass

    def get_sample_idx(self, sample_id: str) -> int:
        """Get the index of a sample in the sample list.

        Args:
            sample_id: Sample ID to look up

        Returns:
            0-based index of the sample

        Raises:
            ValueError: If sample ID not found in the file
        """
        try:
            return self.samples.index(sample_id)
        except ValueError as e:
            raise ValueError(
                f"Sample '{sample_id}' not found. Available samples: {', '.join(self.samples)}"
            ) from e

    def validate(self) -> bool:
        """Validate that the variant source can be accessed.

        Performs basic validation by attempting to access variant_count
        and sample_count properties. Subclasses can override for more
        thorough validation.

        Returns:
            True if validation successful, False otherwise
        """
        try:
            # Try to access basic properties
            _ = self.variant_count
            _ = self.sample_count
            return True
        except Exception:
            return False

    def close(self):  # noqa: B027
        """Close the variant source and release resources.

        Default implementation does nothing. Subclasses should override
        if they need to clean up resources (close file handles, etc.).
        """
        pass

    def __enter__(self) -> "VariantSource":
        """Enter context manager.

        Returns:
            self for use in with statements
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager and clean up resources.

        Args:
            exc_type: Exception type if an error occurred
            exc_val: Exception value if an error occurred
            exc_tb: Exception traceback if an error occurred

        Returns:
            None (does not suppress exceptions)
        """
        self.close()
        return None
