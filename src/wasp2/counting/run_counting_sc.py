import sys
import timeit
import re
import functools
import tempfile
import warnings
from pathlib import Path
from typing import Any, Callable, List, Optional, Union, Iterable, Tuple

import anndata as ad

# local imports
from wasp2.counting.filter_variant_data import (
    vcf_to_bed,
    intersect_vcf_region,
    parse_intersect_region_new,
)
from wasp2.counting.run_counting import tempdir_decorator
from wasp2.counting.count_alleles_sc import make_count_matrix


class WaspCountSC:
    """
    Class to store and process input file paths and parameters for single-cell variant counting.

    This class holds user-specified file paths (BAM, VCF, barcode, and feature files) along with
    additional parameters for variant counting in single-cell experiments. It determines file prefixes,
    creates paths for intermediary files (such as filtered VCF in BED format and intersection files),
    and sets flags for processing gene files versus region files.

    Attributes
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    barcode_file : str
        Path to the barcode file (with one barcode per line).
    feature_file : str
        Path to the feature file (e.g., BED, GTF, or GFF3) used for intersecting variants.
    samples : Optional[List[str]]
        Sample information, converted to a list.
    use_region_names : bool
        Flag indicating whether to use region (or feature) names instead of coordinates.
    out_file : str
        Output file path for the resulting AnnData object (H5AD format).
    temp_loc : str
        Directory for storing intermediary files.
    feature_type : Optional[str]
        Indicates the type of feature file ("regions" or "genes").
    is_gene_file : bool
        True if the feature file is a gene file (GTF/GFF3), else False.
    vcf_prefix : str
        Prefix derived from the VCF file name.
    vcf_bed : str
        Path to the filtered VCF output in BED format.
    intersect_file : str
        Path where the intersection output will be written.
    gtf_bed : Optional[str]
        For gene files: path to the BED file converted from the gene file.
    use_feature_names : bool
        Flag to use feature attributes as region names.
    """

    def __init__(
        self,
        bam_file: str,
        vcf_file: str,
        barcode_file: str,
        feature_file: str,
        samples: Optional[Union[str, List[str]]] = None,
        use_region_names: bool = False,
        out_file: Optional[str] = None,
        temp_loc: Optional[str] = None,
    ) -> None:
        """
        Initialize a WaspCountSC instance with the provided file paths and parameters.

        Parameters
        ----------
        bam_file : str
            Path to the BAM file.
        vcf_file : str
            Path to the VCF file.
        barcode_file : str
            Path to the barcode file (one barcode per line). May be optional.
        feature_file : str
            Path to the feature file (BED, GTF, or GFF3) for intersecting variants.
        samples : str or list of str, optional
            Sample information; ideally a single sample for single-cell experiments.
        use_region_names : bool, optional
            Flag to use region (or feature) names instead of coordinates. Defaults to False.
        out_file : str, optional
            Output file for the resulting AnnData object. Defaults to "allele_counts.h5ad" in the current directory.
        temp_loc : str, optional
            Directory for intermediary files. Defaults to the current working directory.
        """
        # TODO: ALSO ACCEPT .h5

        # User input files
        self.bam_file: str = bam_file
        self.vcf_file: str = vcf_file
        self.barcode_file: str = barcode_file  # Maybe could be optional?
        self.feature_file: str = feature_file
        self.samples: Optional[Union[str, List[str]]] = samples
        self.use_region_names: bool = use_region_names
        self.out_file: Optional[str] = out_file
        self.temp_loc: Optional[str] = temp_loc

        # Optional inputs and outputs?
        # output_sparse_mtx = None
        # SNP OUTPUT?!?!?

        # Make sure samples are turned into a str list (ideally a single sample for single cell)
        if isinstance(self.samples, str):
            # Check if sample file or comma-delimited string
            if Path(self.samples).is_file():
                with open(self.samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]
            else:
                self.samples = [s.strip() for s in self.samples.split(",")]

        # Parse output if not provided
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "allele_counts.h5ad")

        # Failsafe if decorator doesn't create temp_loc
        if self.temp_loc is None:
            self.temp_loc = str(Path.cwd())

        # Parse VCF and intersect files
        vcf_prefix = re.split(r'.vcf|.bcf', Path(self.vcf_file).name)[0]
        self.vcf_prefix: str = vcf_prefix

        # Filtered VCF output
        self.vcf_bed: str = str(Path(self.temp_loc) / f"{vcf_prefix}.bed")

        # Parse feature file
        self.feature_type: Optional[str] = None  # maybe use a boolean flag instead

        if self.feature_file is not None:
            f_ext = "".join(Path(self.feature_file).suffixes)

            if re.search(r'\.(.*Peak|bed)(?:\.gz)?$', f_ext, re.I):
                self.feature_type = "regions"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_regions.bed")
                self.is_gene_file = False
            elif re.search(r'\.g[tf]f(?:\.gz)?$', f_ext, re.I):
                self.feature_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.g[tf]f', Path(self.feature_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_feature_names = True  # Use feature attributes as region names
            elif re.search(r'\.gff3(?:\.gz)?$', f_ext, re.I):
                self.feature_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.gff3', Path(self.feature_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_feature_names = True  # Use feature attributes as feature names
            else:
                self.feature_file = None
                print("invalid ftype")  # Make this raise an error later
        else:
            self.intersect_file = self.vcf_bed

        # TODO UPDATE THIS WHEN I ADD AUTOPARSERS
        if self.is_gene_file:
            # Possible edge case of VCF and GTF prefix conflict
            if self.vcf_bed == self.gtf_bed:
                self.gtf_bed = str(Path(self.temp_loc) / "genes.bed")


@tempdir_decorator
def run_count_variants_sc(
    bam_file: str,
    vcf_file: str,
    barcode_file: str,
    feature_file: Optional[str] = None,
    samples: Optional[Union[str, List[str]]] = None,
    use_region_names: bool = False,
    out_file: Optional[str] = None,
    temp_loc: Optional[str] = None,
) -> None:
    """
    Run the single-cell variant counting pipeline.

    This function processes single-cell variant data by filtering a VCF file (using bcftools),
    intersecting it with a feature file (BED, GTF, or GFF3), and then counting alleles from a BAM file.
    It creates intermediary files and generates an AnnData object containing allele counts.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    barcode_file : str
        Path to the barcode file (one barcode per line).
    feature_file : str, optional
        Path to the feature file (BED, GTF, or GFF3) used for intersecting variants.
    samples : str or list of str, optional
        Sample information; ideally a single sample for single-cell experiments.
    use_region_names : bool, optional
        Flag to use region names (or feature attributes) instead of coordinates.
    out_file : str, optional
        Output file for the resulting AnnData object in H5AD format.
    temp_loc : str, optional
        Directory for intermediary files. If not provided, a temporary directory is used.
    
    Returns
    -------
    None
    """
    # Stores file names in the data class
    count_files = WaspCountSC(
        bam_file,
        vcf_file,
        barcode_file=barcode_file,
        feature_file=feature_file,
        samples=samples,
        use_region_names=use_region_names,
        out_file=out_file,
        temp_loc=temp_loc,
    )
    
    print(*vars(count_files).items(), sep="\n")  # For debugging
    
    # Create intermediary files
    # Maybe change include_gt based on preparse?
    vcf_to_bed(
        vcf_file=count_files.vcf_file,
        out_bed=count_files.vcf_bed,
        samples=count_files.samples,
        include_gt=True,
    )

    intersect_vcf_region(
        vcf_file=count_files.vcf_bed,
        region_file=count_files.feature_file,
        out_file=count_files.intersect_file,
    )

    # TODO: handle use_region_names better
    df = parse_intersect_region_new(
        intersect_file=count_files.intersect_file,
        samples=count_files.samples,
        use_region_names=use_region_names,
        region_col=None,
    )
    
    # TODO: handle case where barcode file contains multiple columns
    with open(count_files.barcode_file, "r") as file:
        bc_dict = {line.rstrip(): i for i, line in enumerate(file)}
    
    # Generate Output AnnData object
    adata = make_count_matrix(
        bam_file=count_files.bam_file,
        df=df,
        bc_dict=bc_dict,
        include_samples=count_files.samples,
    )
    
    # Write outputs
    adata.write_h5ad(count_files.out_file)
    # TODO: include output options, (ie MTX, dense?)
