import sys
import timeit
import re
import functools
import tempfile
import warnings
from pathlib import Path
from typing import Any, Callable, List, Optional, Union

# local imports
from wasp2.counting.filter_variant_data import (
    vcf_to_bed,
    intersect_vcf_region,
    parse_intersect_region,
    parse_intersect_region_new,
)
from wasp2.counting.parse_gene_data import (
    make_gene_data,
    parse_intersect_genes,
    parse_intersect_genes_new,
)
from wasp2.counting.count_alleles import make_count_df


class WaspCountFiles:
    """
    Class to hold and process input file paths and parameters for variant counting.

    This class stores the paths to the BAM and VCF files, region files (if provided),
    sample information, and various output paths. It also determines the type of region file
    (regions vs. genes) based on file extensions and prepares file prefixes for downstream
    processing. For gene files, it sets additional flags and output paths accordingly.

    Attributes
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    region_file : Optional[str]
        Path to the region file (BED, GTF, GFF3), if provided.
    samples : Optional[List[str]]
        List of sample names.
    use_region_names : bool
        Flag to indicate whether to use region names (from features) instead of coordinates.
    out_file : str
        Output file path for the final counts. Defaults to "counts.tsv" in the current working directory.
    temp_loc : str
        Directory for storing intermediary files. Defaults to the current working directory.
    is_gene_file : bool
        Flag set to True if the region file is a gene file (GTF/GFF3).
    gtf_bed : Optional[str]
        Output path for the BED file converted from a gene file.
    vcf_prefix : str
        Prefix derived from the VCF filename.
    vcf_bed : str
        Path for the filtered VCF output in BED format.
    region_type : Optional[str]
        Type of region file ("regions" or "genes").
    intersect_file : str
        Path where the intersection output will be written.
    """

    def __init__(
        self,
        bam_file: str,
        vcf_file: str,
        region_file: Optional[str] = None,
        samples: Optional[Union[str, List[str]]] = None,
        use_region_names: bool = False,
        out_file: Optional[str] = None,
        temp_loc: Optional[str] = None,
    ) -> None:
        """
        Initialize a WaspCountFiles instance with the given file paths and parameters.

        Parameters
        ----------
        bam_file : str
            Path to the BAM file.
        vcf_file : str
            Path to the VCF file.
        region_file : str, optional
            Path to the region file (BED, GTF, or GFF3). Defaults to None.
        samples : str or list of str, optional
            Sample information, either as a file path or a comma-delimited string.
        use_region_names : bool, optional
            Whether to use region names from the region file instead of coordinates. Defaults to False.
        out_file : str, optional
            Output file for counts. Defaults to "counts.tsv" in the current working directory.
        temp_loc : str, optional
            Directory for intermediary files. Defaults to the current working directory.
        """
        # User input files
        self.bam_file: str = bam_file
        self.vcf_file: str = vcf_file
        self.region_file: Optional[str] = region_file
        self.samples: Optional[List[str]] = samples  # May be converted to list later
        self.use_region_names: bool = use_region_names
        self.out_file: Optional[str] = out_file
        self.temp_loc: Optional[str] = temp_loc

        # gtf and gff specific
        self.is_gene_file: bool = False  # check if using gff3/gtf
        self.gtf_bed: Optional[str] = None

        # Make sure samples turned into str list
        if isinstance(self.samples, str):
            # Check if sample file or comma delim string
            if Path(self.samples).is_file():
                with open(self.samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]
            else:
                self.samples = [s.strip() for s in self.samples.split(",")]

        # Parse output?
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "counts.tsv")

        # Failsafe if decorator doesn't create temp_loc
        if self.temp_loc is None:
            self.temp_loc = str(Path.cwd())

        # Parse vcf and intersect
        vcf_prefix = re.split(r'.vcf|.bcf', Path(self.vcf_file).name)[0]
        self.vcf_prefix: str = vcf_prefix

        # Filtered vcf output
        self.vcf_bed: str = str(Path(self.temp_loc) / f"{vcf_prefix}.bed")

        # Parse region file
        self.region_type: Optional[str] = None  # maybe use a boolean flag instead

        if self.region_file is not None:
            f_ext = "".join(Path(self.region_file).suffixes)

            if re.search(r'\.(.*Peak|bed)(?:\.gz)?$', f_ext, re.I):
                self.region_type = "regions"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_regions.bed")
                self.is_gene_file = False
            elif re.search(r'\.g[tf]f(?:\.gz)?$', f_ext, re.I):
                self.region_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.g[tf]f', Path(self.region_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_region_names = True  # Use feature attributes as region names
            elif re.search(r'\.gff3(?:\.gz)?$', f_ext, re.I):
                self.region_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.gff3', Path(self.region_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_region_names = True  # Use feature attributes as region names
            else:
                self.region_file = None
                print("invalid ftype")  # Make this raise an error later

        else:
            self.intersect_file = self.vcf_bed

        # TODO UPDATE THIS WHEN I ADD AUTOPARSERS
        if self.is_gene_file:
            # Possible edge case of vcf and gtf prefix conflict
            if self.vcf_bed == self.gtf_bed:
                self.gtf_bed = str(Path(self.temp_loc) / "genes.bed")


def tempdir_decorator(func: Callable) -> Callable:
    """
    Decorator to ensure a temporary directory is provided for the 'temp_loc' keyword argument.

    If 'temp_loc' is not provided in the keyword arguments when calling the decorated function,
    this decorator creates a temporary directory, passes it as 'temp_loc', and cleans it up afterwards.

    Parameters
    ----------
    func : callable
        The function to decorate.

    Returns
    -------
    callable
        The wrapped function with temporary directory management.
    """
    @functools.wraps(func)
    def tempdir_wrapper(*args, **kwargs):
        if kwargs.get("temp_loc", None) is not None:
            return func(*args, **kwargs)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                kwargs["temp_loc"] = tmpdir
                return func(*args, **kwargs)
    return tempdir_wrapper


@tempdir_decorator
def run_count_variants(
    bam_file: str,
    vcf_file: str,
    region_file: Optional[str] = None,
    samples: Optional[Union[str, List[str]]] = None,
    use_region_names: Optional[bool] = None,
    out_file: Optional[str] = None,
    temp_loc: Optional[str] = None,
    gene_feature: Optional[str] = None,
    gene_attribute: Optional[str] = None,
    gene_parent: Optional[str] = None,
) -> None:
    """
    Run the variant counting pipeline by filtering and intersecting VCF data.

    This function processes the input BAM and VCF files along with optional region and gene-related
    parameters to perform variant counting. It instantiates a WaspCountFiles object, creates intermediary
    files using bcftools and bedtools, parses intersections, and finally counts alleles from the BAM file,
    writing the resulting counts to the specified output file.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    region_file : str, optional
        Path to the region file (e.g., BED, GTF, GFF3) to use for intersecting variants.
    samples : list of str or str, optional
        Sample(s) to filter the VCF file. Can be a comma-delimited string or file path.
    use_region_names : bool, optional
        Flag to use region names (from feature attributes) instead of coordinates.
    out_file : str, optional
        Path where the final counts file will be written.
    temp_loc : str, optional
        Directory for intermediary files. If not provided, a temporary directory is used.
    gene_feature : str, optional
        Feature type (e.g., "exon") for gene counting.
    gene_attribute : str, optional
        Attribute to extract from gene data.
    gene_parent : str, optional
        Parent attribute for grouping in gene data.

    Returns
    -------
    None
    """
    # call the data class
    count_files = WaspCountFiles(
        bam_file=bam_file,
        vcf_file=vcf_file,
        region_file=region_file,
        samples=samples,
        use_region_names=use_region_names,
        out_file=out_file,
        temp_loc=temp_loc,
    )
    
    # print(*vars(count_files).items(), sep="\n") # For debugging
    with_gt = False
    if (count_files.samples is not None) and (len(count_files.samples) == 1):
        with_gt = True
        # temporarily disable for ASE
        # if not count_files.is_gene_file:
        #     with_gt = True
            
    # Create Intermediary Files
    vcf_to_bed(
        vcf_file=count_files.vcf_file,
        out_bed=count_files.vcf_bed,
        samples=count_files.samples,
        include_gt=with_gt,
    )

    # TODO PARSE GENE FEATURES AND ATTRIBUTES
    region_col_name: Optional[str] = None  # Defaults to 'region' as region name
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
    count_df = make_count_df(bam_file=count_files.bam_file, df=df)
    
    # Write counts
    count_df.write_csv(count_files.out_file, has_header=True, separator="\t")
    
    # Should i return for use in analysis pipeline?
    # return count_df
