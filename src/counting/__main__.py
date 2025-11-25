from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from .run_counting import run_count_variants
from .run_counting_sc import run_count_variants_sc

# app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
app = typer.Typer(pretty_exceptions_short=False)

# TODO GOTTA TEST THIS

@app.command()
def count_variants(
    bam: Annotated[str, typer.Argument(help="BAM file")],
    variants: Annotated[str, typer.Argument(help="Variant file (VCF, VCF.GZ, BCF, or PGEN)")],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "-s",
            help=(
                "One or more samples to use in variant file. "
                "Accepts comma delimited string "
                "or file with one sample per line"
            )
            )] = None,
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
                     )] = None,
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
                ),
            )] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files. "
                "Defaults to removing intermediary files using temp directory")
            )] = None,
    use_region_names: Annotated[
        bool,
        typer.Option("--use_region_names",
                     help=(
                         "Use region names instead of coordinates. "
                         "Names are denoted in fourth column of BED. "
                         "Ignored if no name column in file. "
                         "Defaults to using coordinates."
                         )
                     )] = False,
    gene_feature: Annotated[
        Optional[str],
        typer.Option(
            "--gene_feature",
            "--feature",
            "--feat",
            help=(
                "Feature type in gtf/gff3 for counting intersecting SNPs. "
                "Defaults to 'exon' for snp counting")
            )] = None,
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
                "Defaults to '<feature>_id' in gtf and 'ID' in gff3")
            )] = None,
    gene_parent: Annotated[
        Optional[str],
        typer.Option(
            "--gene_parent",
            "--parent",
            "--parent_feature",
            "--parent_attribute",
            help=(
                "Parent attribute in gtf/gff3 for feature used in counting"
                "Defaults to 'transcript_id' in gtf and 'Parent' in gff3")
            )] = None,
    use_rust: Annotated[
        bool,
        typer.Option(
            "--use-rust/--no-rust",
            help=(
                "Use Rust acceleration for BAM counting (requires wasp2_rust extension). "
                "Defaults to True if extension is available.")
            )] = True,
    vcf_bed: Annotated[
        Optional[str],
        typer.Option(
            "--vcf-bed",
            help="Optional precomputed VCF bed file to skip vcf_to_bed."
        )
    ] = None,
    intersect_bed: Annotated[
        Optional[str],
        typer.Option(
            "--intersect-bed",
            help="Optional precomputed intersect bed file to skip bedtools intersect."
        )
    ] = None,

) -> None:

    # Parse sample string
    # print(samples)
    sample_str: Optional[str]
    if samples is not None and len(samples) > 0:
        sample_str = samples[0]
    else:
        sample_str = None

    # print(sample_str)

    # run
    run_count_variants(bam_file=bam,
                       variant_file=variants,
                       region_file=region_file,
                       samples=sample_str,
                       use_region_names=use_region_names,
                       out_file=out_file,
                       temp_loc=temp_loc,
                       gene_feature=gene_feature,
                       gene_attribute=gene_attribute,
                       gene_parent=gene_parent,
                       use_rust=use_rust,
                       precomputed_vcf_bed=vcf_bed,
                       precomputed_intersect=intersect_bed
                       )
    
    # TODO TEST CASES FOR TYPER
    # TODO UNIT TEST NEXT


@app.command()
def count_variants_sc(
    bam: Annotated[str, typer.Argument(help="BAM file")],
    variants: Annotated[str, typer.Argument(help="Variant file (VCF, VCF.GZ, BCF, or PGEN)")],
    barcodes: Annotated[str, typer.Argument(
        help="File with one barcode per line. Used as index")],
    samples: Annotated[
        Optional[List[str]],
        typer.Option(
            "--samples",
            "--sample",
            "--samps",
            "-s",
            help=(
                "One or more samples to use in variant file. "
                "Accepts comma delimited string "
                "or file with one sample per line. "
                "RECOMMENDED TO USE ONE SAMPLE AT A TIME."
            )
            )] = None,
    feature_file: Annotated[
        Optional[str],
        typer.Option("--feature",
                     "--features",
                     "--feat",
                     "-f",
                     "--region",
                     "--regions",
                     "-r",
                     help=("Features used in single-cell experiment. "
                         "Only use variants overlapping features in file. "
                         "Accepts BED or MACS2 formatted .(narrow/broad)Peak files. "
                         "TODO: Implement genes gtf/gff format"
                     )
                     )] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--out",
            "-o",
            help=("Output file to write Anndata allele counts. "
                "H5ad file format. "
                "Defaults to allele_counts.h5ad"
                ),
            )] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files. "
                "Defaults to removing intermediary files using temp directory")
            )] = None
) -> None:

    # Parse sample string
    sample_str: Optional[str]
    if samples is not None and len(samples) > 0:
        sample_str = samples[0]
    else:
        sample_str = None

    # run
    run_count_variants_sc(bam_file=bam,
                          variant_file=variants,
                          barcode_file=barcodes,
                          feature_file=feature_file,
                          samples=sample_str,
                          out_file=out_file,
                          temp_loc=temp_loc
                          )


if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()
