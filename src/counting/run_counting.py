from __future__ import annotations

import functools
import re
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import ParamSpec, TypeVar

from .count_alleles import make_count_df

# local imports
from .filter_variant_data import intersect_vcf_region, parse_intersect_region_new, vcf_to_bed
from .parse_gene_data import make_gene_data, parse_intersect_genes_new

P = ParamSpec("P")
T = TypeVar("T")


class WaspCountFiles:
    """Container for WASP counting pipeline file paths and configuration.

    Manages input/output file paths and parsing logic for the variant counting pipeline.

    Attributes:
        bam_file: Path to the BAM alignment file.
        variant_file: Path to the variant file (VCF, BCF, or PGEN).
        region_file: Optional path to a region file (BED, GTF, or GFF3).
        samples: List of sample IDs to process, or None for all samples.
        use_region_names: Whether to use region names from the region file.
        out_file: Output file path for count results.
        temp_loc: Directory for temporary files.
        is_gene_file: Whether the region file is a gene annotation file.
        gtf_bed: Path to converted GTF/GFF3 BED file, if applicable.
        variant_prefix: Prefix extracted from variant filename.
        vcf_bed: Path to variant BED file.
        skip_vcf_to_bed: Whether to skip VCF-to-BED conversion.
        region_type: Type of regions ('regions' or 'genes').
        intersect_file: Path to intersected variant-region file.
        skip_intersect: Whether to skip intersection step.
    """

    # Class attribute type hints
    bam_file: str
    variant_file: str
    region_file: str | None
    samples: list[str] | None
    use_region_names: bool
    out_file: str
    temp_loc: str
    is_gene_file: bool
    gtf_bed: str | None
    variant_prefix: str
    vcf_bed: str
    skip_vcf_to_bed: bool
    region_type: str | None
    intersect_file: str
    skip_intersect: bool

    def __init__(
        self,
        bam_file: str,
        variant_file: str,
        region_file: str | None = None,
        samples: str | list[str] | None = None,
        use_region_names: bool = False,
        out_file: str | None = None,
        temp_loc: str | None = None,
        precomputed_vcf_bed: str | None = None,
        precomputed_intersect: str | None = None,
    ) -> None:
        # User input files
        self.bam_file = bam_file
        self.variant_file = variant_file
        self.region_file = region_file
        self.use_region_names = use_region_names

        # gtf and gff specific
        self.is_gene_file = False  # check if using gff3/gtf
        self.gtf_bed = None

        # Make sure samples turned into str list
        if isinstance(samples, str):
            # Check if sample file or comma delim string
            if Path(samples).is_file():
                with open(samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]
            else:
                self.samples = [s.strip() for s in samples.split(",")]
        else:
            self.samples = samples

        # parse output?
        self.out_file: str = out_file if out_file is not None else str(Path.cwd() / "counts.tsv")

        # Failsafe if decorator doesnt create temp_loc
        self.temp_loc: str = temp_loc if temp_loc is not None else str(Path.cwd())

        # Parse variant file prefix (handle VCF, BCF, PGEN)
        variant_name = Path(self.variant_file).name
        if variant_name.endswith(".vcf.gz"):
            variant_prefix = variant_name[:-7]  # Remove .vcf.gz
        elif variant_name.endswith(".pgen"):
            variant_prefix = variant_name[:-5]  # Remove .pgen
        else:
            variant_prefix = re.split(r"\.vcf|\.bcf", variant_name)[0]
        self.variant_prefix = variant_prefix

        # Filtered variant output (or precomputed)
        self.vcf_bed = (
            precomputed_vcf_bed
            if precomputed_vcf_bed is not None
            else str(Path(self.temp_loc) / f"{variant_prefix}.bed")
        )
        self.skip_vcf_to_bed = precomputed_vcf_bed is not None

        # Parse region file
        self.region_type = None  # maybe use a boolean flag instead

        if self.region_file is not None:
            f_ext = "".join(Path(self.region_file).suffixes)

            if re.search(r"\.(.*Peak|bed)(?:\.gz)?$", f_ext, re.I):
                self.region_type = "regions"
                self.intersect_file = (
                    precomputed_intersect
                    if precomputed_intersect is not None
                    else str(Path(self.temp_loc) / f"{variant_prefix}_intersect_regions.bed")
                )
                self.is_gene_file = False
            elif re.search(r"\.g[tf]f(?:\.gz)?$", f_ext, re.I):
                self.region_type = "genes"
                self.intersect_file = (
                    precomputed_intersect
                    if precomputed_intersect is not None
                    else str(Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed")
                )
                self.is_gene_file = True
                gtf_prefix = re.split(r".g[tf]f", Path(self.region_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_region_names = True  # Use feature attributes as region names
            elif re.search(r"\.gff3(?:\.gz)?$", f_ext, re.I):
                self.region_type = "genes"
                self.intersect_file = (
                    precomputed_intersect
                    if precomputed_intersect is not None
                    else str(Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed")
                )
                self.is_gene_file = True
                gtf_prefix = re.split(r".gff3", Path(self.region_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_region_names = True  # Use feature attributes as region names
            else:
                raise ValueError(
                    f"Invalid region file type. Expected .bed, .gtf, or .gff3, got: {self.region_file}"
                )

        else:
            # No region file: intersect file defaults to vcf_bed (or provided precomputed)
            self.intersect_file = (
                precomputed_intersect if precomputed_intersect is not None else self.vcf_bed
            )
        self.skip_intersect = precomputed_intersect is not None

        # TODO UPDATE THIS WHEN I ADD AUTOPARSERS
        if self.is_gene_file:
            # Possible edge case of vcf and gtf prefix conflict
            if self.vcf_bed == self.gtf_bed:
                self.gtf_bed = str(Path(self.temp_loc) / "genes.bed")


def tempdir_decorator(func: Callable[P, T]) -> Callable[P, T]:
    """Decorator that creates a temporary directory for the wrapped function.

    If 'temp_loc' is not provided in kwargs, creates a temporary directory
    and passes it to the function. The directory is cleaned up after execution.

    Args:
        func: The function to wrap.

    Returns:
        Wrapped function with automatic temporary directory management.
    """

    @functools.wraps(func)
    def tempdir_wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        if kwargs.get("temp_loc") is not None:
            return func(*args, **kwargs)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                kwargs["temp_loc"] = tmpdir
                return func(*args, **kwargs)

    return tempdir_wrapper


@tempdir_decorator
def run_count_variants(
    bam_file: str,
    variant_file: str,
    region_file: str | None = None,
    samples: str | list[str] | None = None,
    use_region_names: bool = False,
    out_file: str | None = None,
    temp_loc: str | None = None,
    gene_feature: str | None = None,
    gene_attribute: str | None = None,
    gene_parent: str | None = None,
    use_rust: bool = True,
    precomputed_vcf_bed: str | None = None,
    precomputed_intersect: str | None = None,
    include_indels: bool = False,
) -> None:
    """Run the WASP variant counting pipeline.

    Counts allele-specific reads at heterozygous variant positions within
    optional genomic regions.

    Args:
        bam_file: Path to the BAM alignment file.
        variant_file: Path to the variant file (VCF, BCF, or PGEN).
        region_file: Optional path to a region file (BED, GTF, or GFF3).
        samples: Sample ID(s) to process. Can be a single ID, comma-separated
            string, path to a file with one sample per line, or list of IDs.
        use_region_names: Whether to use region names from the region file.
        out_file: Output file path. Defaults to 'counts.tsv' in current directory.
        temp_loc: Directory for temporary files. Auto-created if not provided.
        gene_feature: GTF/GFF3 feature type to extract (e.g., 'gene', 'exon').
        gene_attribute: GTF/GFF3 attribute for region names (e.g., 'gene_name').
        gene_parent: GTF/GFF3 parent attribute for hierarchical features.
        use_rust: Whether to use the Rust backend for counting (faster).
        precomputed_vcf_bed: Path to pre-computed variant BED file (skips conversion).
        precomputed_intersect: Path to pre-computed intersection file.
        include_indels: Whether to include indels in variant counting.

    Returns:
        None. Results are written to out_file.
    """
    # call the data class
    count_files = WaspCountFiles(
        bam_file,
        variant_file,
        region_file=region_file,
        samples=samples,
        use_region_names=use_region_names,
        out_file=out_file,
        temp_loc=temp_loc,
        precomputed_vcf_bed=precomputed_vcf_bed,
        precomputed_intersect=precomputed_intersect,
    )

    # print(*vars(count_files).items(), sep="\n") # For debugging
    with_gt = False
    if (count_files.samples is not None) and (len(count_files.samples) == 1):
        with_gt = True

        # temporarily disable for ASE
        # if not count_files.is_gene_file:
        #     with_gt = True

    # Create Intermediary Files
    if not count_files.skip_vcf_to_bed:
        vcf_to_bed(
            vcf_file=count_files.variant_file,
            out_bed=count_files.vcf_bed,
            samples=count_files.samples,
            include_gt=with_gt,
            include_indels=include_indels,
        )

    # TODO PARSE GENE FEATURES AND ATTRIBUTES
    region_col_name = None  # Defaults to 'region' as region name
    intersect_genes = False

    # region_files is valid to perform intersects
    if count_files.region_file is not None:
        # Check if we need to prepare genes for intersection
        if count_files.gtf_bed is not None:
            # TODO UPDATE THIS WHEN I ADD AUTOPARSERS AND VALIDATORS
            gene_data = make_gene_data(
                gene_file=count_files.region_file,
                out_bed=count_files.gtf_bed,
                feature=gene_feature,
                attribute=gene_attribute,
                parent_attribute=gene_parent,
            )

            regions_to_intersect = count_files.gtf_bed
            region_col_name = gene_data.feature
            intersect_genes = True
        else:
            regions_to_intersect = count_files.region_file
            region_col_name = None  # Defaults to 'region' as region name

        if not count_files.skip_intersect:
            intersect_vcf_region(
                vcf_file=count_files.vcf_bed,
                region_file=regions_to_intersect,
                out_file=count_files.intersect_file,
            )

    # Create Variant Dataframe
    # TODO validate
    if intersect_genes:
        df = parse_intersect_genes_new(
            intersect_file=count_files.intersect_file,
            attribute=gene_data.attribute,
            parent_attribute=gene_data.parent_attribute,
        )
    elif with_gt:
        df = parse_intersect_region_new(
            intersect_file=count_files.intersect_file,
            samples=["GT"],
            use_region_names=count_files.use_region_names,
            region_col=region_col_name,
        )
    else:
        df = parse_intersect_region_new(
            intersect_file=count_files.intersect_file,
            samples=None,
            use_region_names=count_files.use_region_names,
            region_col=region_col_name,
        )
        # df = parse_intersect_region(
        #     intersect_file=count_files.intersect_file,
        #     use_region_names=count_files.use_region_names,
        #     region_col=region_col_name)

    # Should I include a filt bam step???

    # Count
    count_df = make_count_df(bam_file=count_files.bam_file, df=df, use_rust=use_rust)

    # Write counts
    count_df.write_csv(count_files.out_file, include_header=True, separator="\t")

    # Should i return for use in analysis pipeline?
    # return count_df
