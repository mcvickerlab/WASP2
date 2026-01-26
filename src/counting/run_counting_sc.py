"""Single-cell variant counting pipeline."""

from __future__ import annotations

import re
from pathlib import Path

from .count_alleles_sc import make_count_matrix

# local imports
from .filter_variant_data import intersect_vcf_region, parse_intersect_region_new, vcf_to_bed
from .run_counting import tempdir_decorator


class WaspCountSC:
    """Container for single-cell WASP counting pipeline configuration.

    Attributes
    ----------
    bam_file : str
        Path to the BAM alignment file.
    variant_file : str
        Path to the variant file (VCF, BCF, or PGEN).
    barcode_file : str
        Path to cell barcode file.
    feature_file : str | None
        Optional path to feature/region file.
    samples : list[str] | None
        List of sample IDs to process.
    out_file : str
        Output file path for AnnData.
    """

    def __init__(
        self,
        bam_file: str,
        variant_file: str,
        barcode_file: str,
        feature_file: str | None,
        samples: str | list[str] | None = None,
        use_region_names: bool = False,
        out_file: str | None = None,
        temp_loc: str | None = None,
    ) -> None:
        # TODO: ALSO ACCEPT .h5

        # User input files
        self.bam_file = bam_file
        self.variant_file = variant_file
        self.barcode_file = barcode_file  # Maybe could be optional?

        self.feature_file = feature_file
        self.samples = samples
        self.use_region_names = use_region_names
        self.out_file = out_file
        self.temp_loc = temp_loc

        # Optional inputs and outputs?
        # output_sparse_mtx = None
        # SNP OUTPUT?!?!?

        # Make sure samples turned into str list
        # Ideally single sample for single cell
        if isinstance(self.samples, str):
            # Check if sample file or comma delim string
            if Path(self.samples).is_file():
                with open(self.samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]

            else:
                self.samples = [s.strip() for s in self.samples.split(",")]

        # parse output?
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "allele_counts.h5ad")

        # Failsafe if decorator doesnt create temp_loc
        if self.temp_loc is None:
            self.temp_loc = str(Path.cwd())

        # Parse variant file prefix (handle VCF, BCF, PGEN)
        variant_name = Path(self.variant_file).name
        if variant_name.endswith(".vcf.gz"):
            variant_prefix = variant_name[:-7]  # Remove .vcf.gz
        elif variant_name.endswith(".pgen"):
            variant_prefix = variant_name[:-5]  # Remove .pgen
        else:
            variant_prefix = re.split(r"\.vcf|\.bcf", variant_name)[0]
        self.variant_prefix = variant_prefix

        # Filtered variant output
        self.vcf_bed = str(Path(self.temp_loc) / f"{variant_prefix}.bed")

        # Parse feature file
        self.feature_type = None  # maybe use a boolean flag instead

        if self.feature_file is not None:
            f_ext = "".join(Path(self.feature_file).suffixes)

            if re.search(r"\.(.*Peak|bed)(?:\.gz)?$", f_ext, re.I):
                self.feature_type = "regions"
                self.intersect_file = str(
                    Path(self.temp_loc) / f"{variant_prefix}_intersect_regions.bed"
                )
                self.is_gene_file = False
            elif re.search(r"\.g[tf]f(?:\.gz)?$", f_ext, re.I):
                self.feature_type = "genes"
                self.intersect_file = str(
                    Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed"
                )
                self.is_gene_file = True
                gtf_prefix = re.split(r".g[tf]f", Path(self.feature_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_feature_names = True  # Use feature attributes as region names
            elif re.search(r"\.gff3(?:\.gz)?$", f_ext, re.I):
                self.feature_type = "genes"
                self.intersect_file = str(
                    Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed"
                )
                self.is_gene_file = True
                gtf_prefix = re.split(r".gff3", Path(self.feature_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_feature_names = True  # Use feature attributes as feature names
            else:
                raise ValueError(
                    f"Invalid feature file type. Expected .bed, .gtf, or .gff3, got: {self.feature_file}"
                )

        else:
            self.intersect_file = self.vcf_bed

        # TODO UPDATE THIS WHEN I ADD AUTOPARSERS
        if self.is_gene_file:
            # Possible edge case of vcf and gtf prefix conflict
            if self.vcf_bed == self.gtf_bed:
                self.gtf_bed = str(Path(self.temp_loc) / "genes.bed")


@tempdir_decorator
def run_count_variants_sc(
    bam_file: str,
    variant_file: str,
    barcode_file: str,
    feature_file: str | None = None,
    samples: str | list[str] | None = None,
    use_region_names: bool = False,
    out_file: str | None = None,
    temp_loc: str | None = None,
) -> None:
    """Run single-cell variant counting pipeline.

    Parameters
    ----------
    bam_file : str
        Path to the BAM alignment file with cell barcodes.
    variant_file : str
        Path to the variant file (VCF, BCF, or PGEN).
    barcode_file : str
        Path to cell barcode file (one barcode per line).
    feature_file : str | None, optional
        Path to feature/region file (BED, GTF, or GFF3).
    samples : str | list[str] | None, optional
        Sample ID(s) to process.
    use_region_names : bool, optional
        Whether to use region names from the feature file.
    out_file : str | None, optional
        Output file path. Defaults to 'allele_counts.h5ad'.
    temp_loc : str | None, optional
        Directory for temporary files.

    Returns
    -------
    None
        Results are written to out_file as AnnData.
    """
    # Stores file names
    count_files = WaspCountSC(
        bam_file,
        variant_file,
        barcode_file=barcode_file,
        feature_file=feature_file,
        samples=samples,
        use_region_names=use_region_names,
        out_file=out_file,
        temp_loc=temp_loc,
    )

    # Create intermediary files
    # Maybe change include_gt based on preparse?
    vcf_to_bed(
        vcf_file=count_files.variant_file,
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
    with open(count_files.barcode_file) as file:
        bc_dict = {line.rstrip(): i for i, line in enumerate(file)}

    # Generate Output
    adata = make_count_matrix(
        bam_file=count_files.bam_file, df=df, bc_dict=bc_dict, include_samples=count_files.samples
    )

    # Write outputs
    adata.write_h5ad(count_files.out_file)
    # TODO: include output options, (ie MTX, dense?)
