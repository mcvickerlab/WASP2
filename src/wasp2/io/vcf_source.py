"""
VCF/BCF reader implementation for WASP2.

This module provides VCFSource, a VariantSource implementation for reading
VCF and BCF files using pysam. Supports both plain and compressed formats,
with optional indexing for region queries.

When available, uses Rust acceleration (wasp2_rust) for VCF → BED conversion
which is 5-6x faster than bcftools subprocess.
"""

import os
import subprocess
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import pysam

from .variant_source import (
    Genotype,
    Variant,
    VariantGenotype,
    VariantSource,
)

# Try to import Rust acceleration
try:
    from wasp2_rust import vcf_to_bed_py as rust_vcf_to_bed
    RUST_VCF_AVAILABLE = True
except ImportError:
    RUST_VCF_AVAILABLE = False


@VariantSource.register('vcf', 'vcf.gz', 'vcf.bgz', 'bcf', 'bcf.gz')
class VCFSource(VariantSource):
    """VariantSource implementation for VCF/BCF files.

    Reads variant data from VCF (Variant Call Format) and BCF (binary VCF) files
    using pysam/htslib. Supports both plain and compressed formats (.vcf, .vcf.gz, .bcf),
    and can leverage tabix/CSI indexes for efficient region queries.

    The class handles:
    - Standard VCF/BCF parsing
    - Genotype extraction and conversion to Genotype enum
    - Sample-specific filtering
    - Heterozygous-only filtering
    - Region queries (if indexed)
    - BED format export using bcftools for efficiency

    Attributes:
        path: Path to the VCF/BCF file
        vcf: pysam.VariantFile handle
        _samples: Cached list of sample IDs
        _variant_count: Cached variant count (lazy computed)

    Example:
        >>> with VCFSource("variants.vcf.gz") as vcf:
        ...     for vg in vcf.iter_variants(het_only=True):
        ...         print(f"{vg.variant.chrom}:{vg.variant.pos}")
    """

    def __init__(self, path: str, **kwargs):
        """Initialize VCF source.

        Args:
            path: Path to VCF/BCF file (str or Path-like)
            **kwargs: Additional arguments (reserved for future use)

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file cannot be opened or parsed
        """
        self.path = Path(path)

        # Open VCF file with pysam
        try:
            self.vcf = pysam.VariantFile(str(self.path))
        except (OSError, ValueError) as e:
            raise ValueError(f"Failed to open VCF file {self.path}: {e}")

        # Cache samples from header
        self._samples = list(self.vcf.header.samples)

        # Lazy-computed variant count
        self._variant_count: Optional[int] = None

    @property
    def samples(self) -> List[str]:
        """Get list of sample IDs from VCF header.

        Returns:
            List of sample ID strings in file order
        """
        return self._samples

    @property
    def variant_count(self) -> int:
        """Get total number of variants in the file.

        Counts variants by iterating through the file. Result is cached
        for subsequent calls.

        Returns:
            Total number of variants
        """
        if self._variant_count is None:
            # Count variants by iterating through file
            count = 0
            for _ in self.vcf.fetch():
                count += 1
            self._variant_count = count

            # Reset iterator for future use
            self.vcf.close()
            self.vcf = pysam.VariantFile(str(self.path))

        return self._variant_count

    @property
    def sample_count(self) -> int:
        """Get total number of samples.

        Returns:
            Total number of samples
        """
        return len(self._samples)

    def iter_variants(
        self,
        samples: Optional[List[str]] = None,
        het_only: bool = False
    ) -> Iterator[VariantGenotype]:
        """Iterate over variants with optional filtering.

        Yields one VariantGenotype per variant for the first sample in the list
        (or first sample in file if samples=None).

        Args:
            samples: Optional list of sample IDs. If None, uses first sample.
                    Currently only supports single sample iteration.
            het_only: If True, only yield heterozygous variants

        Yields:
            VariantGenotype objects for each variant

        Example:
            >>> for vg in source.iter_variants(samples=["sample1"], het_only=True):
            ...     print(vg.variant.pos, vg.genotype)
        """
        # Determine which sample to iterate
        if samples is None:
            target_samples = [self._samples[0]] if self._samples else []
        else:
            # Validate samples exist
            for s in samples:
                if s not in self._samples:
                    raise ValueError(f"Sample '{s}' not found in VCF")
            target_samples = samples

        if not target_samples:
            return

        # Currently support single sample iteration
        # (multi-sample would yield multiple VariantGenotype per variant)
        sample_id = target_samples[0]

        # Iterate through VCF records
        for record in self.vcf.fetch():
            # Get sample genotype
            sample_data = record.samples[sample_id]
            gt = sample_data.get('GT', None)

            if gt is None or None in gt:
                # Missing genotype
                genotype = Genotype.MISSING
            else:
                # Parse GT tuple
                genotype = self._parse_gt(gt)

            # Filter by het_only if requested
            if het_only and genotype != Genotype.HET:
                continue

            # Create Variant object (use first ALT if multi-allelic)
            alt = record.alts[0] if record.alts else record.ref
            variant = Variant(
                chrom=record.chrom,
                pos=record.pos,
                ref=record.ref,
                alt=alt,
                id=record.id
            )

            # Get allele sequences
            allele1, allele2 = self._get_alleles(record, gt)

            yield VariantGenotype(
                variant=variant,
                genotype=genotype,
                allele1=allele1,
                allele2=allele2
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
        # Validate sample exists
        if sample not in self._samples:
            raise ValueError(f"Sample '{sample}' not found in VCF")

        # Query the position
        try:
            records = list(self.vcf.fetch(chrom, pos - 1, pos))
        except (OSError, ValueError) as e:
            raise ValueError(f"Failed to query position {chrom}:{pos}: {e}")

        if not records:
            raise ValueError(f"No variant found at {chrom}:{pos}")

        # Get genotype from first matching record
        record = records[0]
        sample_data = record.samples[sample]
        gt = sample_data.get('GT', None)

        if gt is None or None in gt:
            return Genotype.MISSING

        return self._parse_gt(gt)

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        samples: Optional[List[str]] = None
    ) -> Iterator[VariantGenotype]:
        """Query variants in a genomic region.

        Requires the VCF to be indexed (.tbi or .csi). Uses 1-based inclusive
        coordinates (VCF standard).

        Args:
            chrom: Chromosome name
            start: 1-based start position (inclusive)
            end: 1-based end position (inclusive)
            samples: Optional list of sample IDs. If None, uses first sample.

        Yields:
            VariantGenotype objects in the region

        Raises:
            ValueError: If the file is not indexed or region is invalid
        """
        # Determine target sample
        if samples is None:
            target_samples = [self._samples[0]] if self._samples else []
        else:
            for s in samples:
                if s not in self._samples:
                    raise ValueError(f"Sample '{s}' not found in VCF")
            target_samples = samples

        if not target_samples:
            return

        sample_id = target_samples[0]

        # Query region (pysam uses 0-based coordinates for fetch)
        try:
            records = self.vcf.fetch(chrom, start - 1, end)
        except (OSError, ValueError) as e:
            raise ValueError(
                f"Failed to query region {chrom}:{start}-{end}. "
                f"File may not be indexed: {e}"
            )

        # Yield VariantGenotype for each record
        for record in records:
            sample_data = record.samples[sample_id]
            gt = sample_data.get('GT', None)

            if gt is None or None in gt:
                genotype = Genotype.MISSING
            else:
                genotype = self._parse_gt(gt)

            # Create Variant (use first ALT)
            alt = record.alts[0] if record.alts else record.ref
            variant = Variant(
                chrom=record.chrom,
                pos=record.pos,
                ref=record.ref,
                alt=alt,
                id=record.id
            )

            allele1, allele2 = self._get_alleles(record, gt)

            yield VariantGenotype(
                variant=variant,
                genotype=genotype,
                allele1=allele1,
                allele2=allele2
            )

    def to_bed(
        self,
        output: Path,
        samples: Optional[List[str]] = None,
        het_only: bool = True,
        include_genotypes: bool = True,
        include_indels: bool = False,
        max_indel_len: int = 10
    ) -> Path:
        """Export variants to BED format file.

        Uses Rust acceleration when available (5-6x faster), falls back to
        bcftools subprocess. BED format uses 0-based start, 1-based end coordinates.

        Format:
        - Without genotypes: chrom\\tstart\\tend\\tref\\talt
        - With genotypes: chrom\\tstart\\tend\\tref\\talt\\tgenotype

        Args:
            output: Output BED file path
            samples: Optional list of sample IDs to include
            het_only: If True, only export heterozygous variants
            include_genotypes: If True, include genotype column(s)
            include_indels: If True, include indels in addition to SNPs
            max_indel_len: Maximum indel length (bp) to include

        Returns:
            Path to the created BED file

        Raises:
            IOError: If conversion fails or file cannot be written
            ValueError: If samples not found
        """
        # Validate samples if provided
        if samples is not None:
            for s in samples:
                if s not in self._samples:
                    raise ValueError(f"Sample '{s}' not found in VCF")

        # Try Rust acceleration first (5-6x faster than bcftools)
        use_rust = (
            RUST_VCF_AVAILABLE
            and os.environ.get("WASP2_DISABLE_RUST") != "1"
        )

        if use_rust:
            try:
                rust_vcf_to_bed(
                    str(self.path),
                    str(output),
                    samples=samples,
                    het_only=het_only,
                    include_indels=include_indels,
                    max_indel_len=max_indel_len,
                    include_genotypes=include_genotypes,
                )
                return output
            except Exception as e:
                print(f"Rust vcf_to_bed failed: {e}, falling back to bcftools")

        # Fallback to bcftools subprocess
        return self._to_bed_bcftools(
            output, samples, het_only, include_genotypes,
            include_indels, max_indel_len
        )

    def _to_bed_bcftools(
        self,
        output: Path,
        samples: Optional[List[str]],
        het_only: bool,
        include_genotypes: bool,
        include_indels: bool,
        max_indel_len: int
    ) -> Path:
        """Export variants to BED using bcftools subprocess (fallback).

        This is the original implementation using bcftools.
        Note: Multi-allelic sites are now included (removed -m2 -M2 filter)
        to match bcftools -g het behavior used by WASP2-Python benchmark.
        """
        # Build bcftools commands based on parameters
        # NOTE: Removed -m2 -M2 biallelic filter to include multi-allelic het sites

        # Base view command
        view_cmd = [
            "bcftools", "view", str(self.path),
        ]

        # Add variant type filter
        if include_indels:
            view_cmd.extend(["-v", "snps,indels"])  # Both SNPs and indels
            # Add indel length filter (max absolute difference in allele lengths)
            # This filters indels where |len(ALT) - len(REF)| > max_indel_len
            view_cmd.extend(["-i", f'strlen(REF)-strlen(ALT)<={max_indel_len} && strlen(ALT)-strlen(REF)<={max_indel_len}'])
        else:
            view_cmd.extend(["-v", "snps"])  # SNPs only (backward compatible)

        view_cmd.append("-Ou")  # uncompressed BCF output

        # Build query command
        query_cmd = [
            "bcftools", "query",
            "-o", str(output),
            "-f"
        ]

        # Configure based on samples and het_only
        if samples is None:
            # No samples: drop genotypes
            view_cmd.append("--drop-genotypes")
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")

            view_process = subprocess.run(
                view_cmd, stdout=subprocess.PIPE, check=True
            )
        else:
            samples_arg = ",".join(samples)
            num_samples = len(samples)

            if num_samples > 1:
                # Multi-sample: filter to variants with at least one non-ref allele
                view_cmd.extend([
                    "-s", samples_arg,
                    "--min-ac", "1",
                    "--max-ac", str((num_samples * 2) - 1)
                ])
                view_process = subprocess.run(
                    view_cmd, stdout=subprocess.PIPE, check=True
                )
            else:
                # Single sample
                view_cmd.extend(["-s", samples_arg])
                subset_process = subprocess.run(
                    view_cmd, stdout=subprocess.PIPE, check=True
                )

                if het_only:
                    # Filter to het genotypes
                    het_view_cmd = ["bcftools", "view", "--genotype", "het", "-Ou"]
                    view_process = subprocess.run(
                        het_view_cmd,
                        input=subset_process.stdout,
                        stdout=subprocess.PIPE,
                        check=True
                    )
                else:
                    view_process = subset_process

            # Add genotype column if requested
            if include_genotypes:
                query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%TGT]\n")
            else:
                query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")

        # Run query command
        try:
            subprocess.run(
                query_cmd,
                input=view_process.stdout,
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise IOError(f"bcftools failed: {e}")

        return output

    def _parse_gt(self, gt_tuple: Tuple[int, ...]) -> Genotype:
        """Convert pysam GT tuple to Genotype enum.

        Args:
            gt_tuple: Genotype tuple from pysam (e.g., (0, 1), (1, 1))

        Returns:
            Genotype enum value

        Examples:
            >>> _parse_gt((0, 0))  # 0/0
            Genotype.HOM_REF
            >>> _parse_gt((0, 1))  # 0/1
            Genotype.HET
            >>> _parse_gt((1, 1))  # 1/1
            Genotype.HOM_ALT
        """
        if None in gt_tuple:
            return Genotype.MISSING

        # Count number of alt alleles
        num_alts = sum(1 for allele in gt_tuple if allele > 0)

        if num_alts == 0:
            return Genotype.HOM_REF
        elif num_alts == len(gt_tuple):
            return Genotype.HOM_ALT
        else:
            return Genotype.HET

    def _get_alleles(
        self, record: pysam.VariantRecord, gt: Optional[Tuple[int, ...]]
    ) -> Tuple[Optional[str], Optional[str]]:
        """Get allele sequences from genotype.

        Args:
            record: pysam VariantRecord
            gt: Genotype tuple (e.g., (0, 1))

        Returns:
            Tuple of (allele1, allele2) sequences

        Examples:
            >>> record.ref = "A"
            >>> record.alts = ["G"]
            >>> _get_alleles(record, (0, 1))
            ("A", "G")
        """
        if gt is None or None in gt:
            return None, None

        alleles = [record.ref] + list(record.alts if record.alts else [])

        try:
            allele1 = alleles[gt[0]] if gt[0] < len(alleles) else None
            allele2 = alleles[gt[1]] if len(gt) > 1 and gt[1] < len(alleles) else None
            return allele1, allele2
        except (IndexError, TypeError):
            return None, None

    def close(self):
        """Close the pysam VariantFile handle.

        Releases file resources. Should be called when done with the source,
        or use context manager protocol.
        """
        if hasattr(self, 'vcf') and self.vcf is not None:
            self.vcf.close()
