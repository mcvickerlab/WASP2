from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from wasp2.counting.run_counting import run_count_variants
from wasp2.counting.run_counting_sc import run_count_variants_sc

# app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
app = typer.Typer(pretty_exceptions_short=False)


@app.command()
def count_variants(
    bam: Annotated[str, typer.Argument(help="Bam File")],
    vcf: Annotated[str, typer.Argument(help="VCF File")],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "--samps",
            "-s",
            help=(
                "One or more samples to use in VCF. "
                "Accepts comma delimited string "
                "or file with one sample per line"
            )
        )
    ] = None,
    region_file: Annotated[
        Optional[str],
        typer.Option("--region",
                     "--regions",
                     "--region_file",
                     "--regions_file",
                     "-r",
                     help=(
                         "Only use variants overlapping regions in file. "
                         "Accepts BED or MACS2 formatted .(narrow/broad)Peak files. "
                     )
        )
    ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--out",
            "-o",
            help=(
                "Output file for counts. "
                "Defaults to counts.tsv"
            )
        )
    ] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files. "
                "Defaults to removing intermediary files using temp directory"
            )
        )
    ] = None,
    use_region_names: Annotated[
        bool,
        typer.Option("--use_region_names",
                     help=(
                         "Use region names instead of coordinates. "
                         "Names are denoted in fourth column of BED. "
                         "Ignored if no name column in file. "
                         "Defaults to using coordinates."
                     )
        )
    ] = False,
    gene_feature: Annotated[
        Optional[str],
        typer.Option(
            "--gene_feature",
            "--feature",
            "--feat",
            help=(
                "Feature type in gtf/gff3 for counting intersecting SNPs. "
                "Defaults to 'exon' for snp counting"
            )
        )
    ] = None,
    gene_attribute: Annotated[
        Optional[str],
        typer.Option(
            "--gene_attribute",
            "--attribute",
            "--attributes",
            "--attrs",
            "--attr",
            help=(
                "Attribute name from gtf/gff3 attribute column to use as ID. "
                "Defaults to '<feature>_id' in gtf and 'ID' in gff3"
            )
        )
    ] = None,
    gene_parent: Annotated[
        Optional[str],
        typer.Option(
            "--gene_parent",
            "--parent",
            "--parent_feature",
            "--parent_attribute",
            help=(
                "Parent attribute in gtf/gff3 for feature used in counting. "
                "Defaults to 'transcript_id' in gtf and 'Parent' in gff3"
            )
        )
    ] = None,
) -> None:
    """
    Count variants by processing a BAM and VCF file, filtering by sample, and intersecting with regions.

    This command runs the variant counting pipeline by calling :func:`run_count_variants` with the specified
    parameters. It supports filtering the VCF file by sample, limiting variants to those overlapping regions,
    and optionally including gene-related parameters for SNP counting.

    Parameters
    ----------
    bam : str
        Path to the BAM file.
    vcf : str
        Path to the VCF file.
    samples : Optional[List[str]]
        One or more samples to use in VCF. Accepts a comma-delimited string or a file with one sample per line.
    region_file : Optional[str]
        Path to a file containing regions (BED or MACS2 formatted peaks) to restrict variant counting.
    out_file : Optional[str]
        Output file for counts. Defaults to "counts.tsv" if not provided.
    temp_loc : Optional[str]
        Directory for keeping intermediary files. Defaults to using the system temporary directory.
    use_region_names : bool
        If True, use region names (from the fourth column of the BED file) instead of coordinates.
    gene_feature : Optional[str]
        Feature type in gtf/gff3 for counting intersecting SNPs. Defaults to 'exon'.
    gene_attribute : Optional[str]
        Attribute name from the gtf/gff3 attribute column to use as ID. Defaults to "<feature>_id" in gtf and "ID" in gff3.
    gene_parent : Optional[str]
        Parent attribute in gtf/gff3 for the feature used in counting. Defaults to "transcript_id" in gtf and "Parent" in gff3.

    Returns
    -------
    None

    Examples
    --------
    >>> count_variants("example.bam", "example.vcf", samples=["Sample1"], region_file="regions.bed")
    """
    # Parse sample string
    # print(samples)
    if samples and len(samples) > 0:
        samples = samples[0]
    else:
        samples = None
    
    # print(samples)
    
    # run
    run_count_variants(
        bam_file=bam,
        vcf_file=vcf,
        region_file=region_file,
        samples=samples,
        use_region_names=use_region_names,
        out_file=out_file,
        temp_loc=temp_loc,
        gene_feature=gene_feature,
        gene_attribute=gene_attribute,
        gene_parent=gene_parent,
    )
    
    # TODO TEST CASES FOR TYPER
    # TODO UNIT TEST NEXT


@app.command()
def count_variants_sc(
    bam: Annotated[str, typer.Argument(help="Bam File")],
    vcf: Annotated[str, typer.Argument(help="VCF File")],
    barcodes: Annotated[str, typer.Argument(help="File with one barcode per line. Used as index")],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "--samps",
            "-s",
            help=(
                "One or more samples to use in VCF. "
                "Accepts comma delimited string "
                "or file with one sample per line. "
                "RECOMMENDED TO USE ONE SAMPLE AT A TIME."
            )
        )
    ] = None,
    feature_file: Annotated[
        Optional[str],
        typer.Option(
            "--feature",
            "--features",
            "--feat",
            "-f",
            "--region",
            "--regions",
            "-r",
            help=(
                "Features used in single-cell experiment. "
                "Only use variants overlapping features in file. "
                "Accepts BED or MACS2 formatted .(narrow/broad)Peak files. "
                "TODO: Implement genes gtf/gff format"
            )
        )
    ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--out",
            "-o",
            help=(
                "Output file to write Anndata allele counts. "
                "H5ad file format. "
                "Defaults to allele_counts.h5ad"
            )
        )
    ] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files. "
                "Defaults to removing intermediary files using temp directory"
            )
        )
    ] = None,
) -> None:
    """
    Count variants for single-cell data and create an AnnData object with allele counts.

    This command runs the single-cell variant counting pipeline by calling
    :func:`run_count_variants_sc` with the provided parameters. It processes the input BAM and VCF files,
    uses the barcode file as an index, and optionally restricts counting to features specified in a BED file.
    
    Parameters
    ----------
    bam : str
        Path to the BAM file.
    vcf : str
        Path to the VCF file.
    barcodes : str
        Path to the file containing one barcode per line.
    samples : Optional[List[str]]
        One or more samples to use in the VCF. Accepts a comma-delimited string or file.
        RECOMMENDED TO USE ONE SAMPLE AT A TIME.
    feature_file : Optional[str]
        Path to a file containing features (BED or MACS2 formatted peaks) used in the single-cell experiment.
        (TODO: Implement support for genes in gtf/gff format.)
    out_file : Optional[str]
        Output file to write the Anndata allele counts. Defaults to "allele_counts.h5ad" if not provided.
    temp_loc : Optional[str]
        Directory for storing intermediary files. Defaults to using a temporary directory.
    
    Returns
    -------
    None

    Examples
    --------
    >>> count_variants_sc("sc.bam", "sc.vcf", "barcodes.txt", samples=["Sample1"], feature_file="features.bed")
    """
    # Parse sample string
    if samples and len(samples) > 0:
        samples = samples[0]
    else:
        samples = None

    # run
    run_count_variants_sc(
        bam_file=bam,
        vcf_file=vcf,
        barcode_file=barcodes,
        feature_file=feature_file,
        samples=samples,
        out_file=out_file,
        temp_loc=temp_loc,
    )


if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()
