"""
Compatibility module for bridging legacy vcf_to_bed with VariantSource.

This module provides backward-compatible functions that can use either:
1. The new VariantSource interface (for VCF, PGEN, etc.)
2. The legacy bcftools subprocess approach (fallback)

The function signatures match the existing vcf_to_bed() in mapping and counting
modules, making it a drop-in replacement.
"""

import subprocess
from pathlib import Path

from .variant_source import VariantSource


def variants_to_bed(
    variant_file: str | Path,
    out_bed: str | Path,
    samples: list[str] | None = None,
    include_gt: bool = True,
    het_only: bool = True,
    use_legacy: bool = False,
    include_indels: bool = False,
    max_indel_len: int = 10,
) -> Path:
    """Convert variant file to BED format.

    This is a unified interface that works with VCF, VCF.GZ, or PGEN files.
    It uses the VariantSource interface when possible, with fallback to
    bcftools for legacy compatibility.

    Args:
        variant_file: Path to variant file (VCF, VCF.GZ, BCF, or PGEN)
        out_bed: Output BED file path
        samples: List of sample IDs to include. If None, no sample filtering.
        include_gt: Include genotype column(s) in output
        het_only: Only include heterozygous sites (when samples specified)
        use_legacy: Force use of legacy bcftools approach (VCF only)
        include_indels: Include indels in addition to SNPs
        max_indel_len: Maximum indel length (bp) to include

    Returns:
        Path to the output BED file

    Note:
        When samples are specified and het_only=True, only heterozygous
        sites for those samples are output.
    """
    variant_file = Path(variant_file)
    out_bed = Path(out_bed)

    # Detect format
    suffix = variant_file.suffix.lower()
    if suffix == ".gz":
        # Check for .vcf.gz
        if variant_file.stem.lower().endswith(".vcf"):
            suffix = ".vcf.gz"
        else:
            suffix = ".gz"

    # Use legacy for VCF when explicitly requested
    if use_legacy and suffix in (".vcf", ".vcf.gz", ".bcf"):
        return _vcf_to_bed_bcftools(
            vcf_file=variant_file,
            out_bed=out_bed,
            samples=samples,
            include_gt=include_gt,
            include_indels=include_indels,
            max_indel_len=max_indel_len,
        )

    # Use VariantSource for all formats
    with VariantSource.open(variant_file) as source:
        source.to_bed(
            out_bed,
            samples=samples,
            het_only=het_only if samples else False,
            include_genotypes=include_gt,
            include_indels=include_indels,
            max_indel_len=max_indel_len,
        )

    return out_bed


def _vcf_to_bed_bcftools(
    vcf_file: str | Path,
    out_bed: str | Path,
    samples: list[str] | None = None,
    include_gt: bool = True,
    include_indels: bool = False,
    max_indel_len: int = 10,
) -> Path:
    """Legacy vcf_to_bed using bcftools subprocess.

    This is the original implementation for backward compatibility.
    Prefer variants_to_bed() which uses VariantSource.

    Note: Multi-allelic sites are now included (removed -m2 -M2 filter)
    to match bcftools -g het behavior used by WASP2-Python benchmark.

    Args:
        vcf_file: Path to VCF/VCF.GZ/BCF file
        out_bed: Output BED file path
        samples: List of sample IDs to filter
        include_gt: Include genotype column in output
        include_indels: Include indels in addition to SNPs
        max_indel_len: Maximum indel length (bp) to include

    Returns:
        Path to output BED file
    """
    vcf_file = Path(vcf_file)
    out_bed = Path(out_bed)

    # Base commands - NOTE: Removed -m2 -M2 to include multi-allelic het sites
    view_cmd = [
        "bcftools",
        "view",
        str(vcf_file),
    ]

    # Add variant type filter
    if include_indels:
        view_cmd.extend(["-v", "snps,indels"])
        # Add indel length filter
        view_cmd.extend(
            [
                "-i",
                f"strlen(REF)-strlen(ALT)<={max_indel_len} && strlen(ALT)-strlen(REF)<={max_indel_len}",
            ]
        )
    else:
        view_cmd.extend(["-v", "snps"])

    view_cmd.append("-Ou")

    query_cmd = ["bcftools", "query", "-o", str(out_bed), "-f"]

    # Parse based on num samples
    if samples is None:
        # No samples - drop genotypes
        view_cmd.append("--drop-genotypes")
        query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")
        view_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
    else:
        # With samples
        samples_arg = ",".join(samples)
        num_samples = len(samples)

        if num_samples > 1:
            # Multi-sample: filter to sites with at least one het
            view_cmd.extend(
                ["-s", samples_arg, "--min-ac", "1", "--max-ac", str((num_samples * 2) - 1)]
            )
            view_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
        else:
            # Single sample: subset then filter to het
            view_cmd.extend(["-s", samples_arg])
            subset_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)

            # Get het genotypes only
            het_cmd = ["bcftools", "view", "--genotype", "het", "-Ou"]
            view_process = subprocess.run(
                het_cmd, input=subset_process.stdout, stdout=subprocess.PIPE, check=True
            )

        # Format string based on include_gt
        if include_gt:
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%GT]\n")
        else:
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")

    # Run query
    subprocess.run(query_cmd, input=view_process.stdout, check=True)

    return out_bed


# Alias for backward compatibility
vcf_to_bed = _vcf_to_bed_bcftools
