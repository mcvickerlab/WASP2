from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from run_analysis import run_ai_analysis
from run_analysis_sc import run_ai_analysis_sc
from run_compare_ai import run_ai_comparison

# app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
app = typer.Typer(pretty_exceptions_short=False)


# TODO GOTTA TEST THIS

# What should i name this?
@app.command()
def find_imbalance(
    counts: Annotated[str, typer.Argument(help="Count File")],
    min: Annotated[
        Optional[int],
        typer.Option(
            "--min",
            "--min_count",
            help=("Minimum allele count for measuring imbalance."
                  " (Default: 10)"
                  )
            )
        ] = None,
    pseudocount: Annotated[
        Optional[int],
        typer.Option(
            "-p",
            "--ps",
            "--pseudo",
            "--pseudocount",
            help=("Pseudocount added when measuring allelic imbalance. "
                  "(Default: 1)"
                  )
            )
        ] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help=(
                "Output file for analysis. "
                "Defaults to ai_results.tsv"
                ),
            )] = None,
    phased: Annotated[
        Optional[bool],
        typer.Option(
            "--phased",
            help=("Calculate allelic imbalance using the phased haplotype model. "
                  "Genotype info must phased and included in allelic count data!"
                  "\nBy default, calculates unphased AI assuming equal liklihood for each haplotype."
                  )
            )] = False,
    model: Annotated[
        Optional[str],
        typer.Option(
            "-m",
            "--model",
            help=(
                "Model used for measuring optimization parameter when finding imbalance. "
                "HIGHLY RECOMMENDED TO LEAVE AS DEFAULT FOR SINGLE DISPERSION MODEL. "
                "Choice of 'single' or 'linear'. "
                "(Default: 'single')"
                ),
            )] = None,
    region_col: Annotated[
        Optional[str],
        typer.Option(
            "--region_col",
            help=(
                "Name of region column for current data..."
                "'region' for ATAC-seq. " 
                "Attribute name for RNA-seq."
                "(Default: Auto-parses if none provided)"
                ),
            )] = None,
    groupby: Annotated[
        Optional[str],
        typer.Option(
            "--groupby",
            "--group",
            "--parent_col",
            "--parent",
            help=(
                "Report allelic imbalance by parent group instead of feature level in RNA-seq counts. \n"
                "Name of parent column. Not valid if no parent column or if using ATAC-seq peaks. \n"
                "(Default: Report by feature level instead of parent level)"
                ),
            )] = None,
    
):
    
    # Run
    run_ai_analysis(count_file=counts,
                    min_count=min,
                    model=model,
                    pseudocount=pseudocount,
                    phased=phased,
                    out_file=out_file,
                    region_col=region_col,
                    groupby=groupby)
    
    # TODO TEST CASES FOR TYPER


@app.command()
def find_imbalance_sc(
    counts: Annotated[str, typer.Argument(help="Count File")],
    bc_map: Annotated[str, typer.Argument(
        help=(
            "Two Column TSV file mapping specific barcodes to some grouping/celltype. "
            "Each line following format [BARCODE]\t[GROUP]"
            )
        )
                      ],
    min: Annotated[
        Optional[int],
        typer.Option(
            "--min",
            "--min_count",
            help=("Minimum allele count per region for measuring imbalance."
                  " (Default: 10)"
                  )
            )
        ] = None,
    pseudocount: Annotated[
        Optional[int],
        typer.Option(
            "-p",
            "--ps",
            "--pseudo",
            "--pseudocount",
            help=("Pseudocount added when measuring allelic imbalance. "
                  "(Default: 1)"
                  )
            )
        ] = None,
    sample: Annotated[
        Optional[str],
        typer.Option(
            "--sample",
            "--samp",
            "-s",
            help=(
                "Use heterozygous genotypes for this sample in count file. "
                "Automatically parses if data contains 0 or 1 sample. "
                "REQUIRED IF COUNT DATA CONTAINS MULTIPLE SAMPLES."
                ),
            )] = None,
    groups: Annotated[
        Optional[List[str]],
        typer.Option(
            "--groups",
            "--group",
            "--celltypes",
            "--g",
            help=(
                "Specific groups in barcode file/bc_map to analyze allelic imbalance in. "
                "Uses all groups in barcode file/bc_map by default."
                ),
            )] = None,
    phased: Annotated[
        Optional[bool],
        typer.Option(
            "--phased/--unphased",
            help=(
                "If genotypes are phased use phasing information to measure imbalance.\n"
                "Otherwise or if --unphased selected, assume all haplotypes are equally likely during analysis.\n"
                "Autoparses genotype data by default if not denoted."
                ))] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help=(
                "Output file for analysis. "
                "Defaults to ai_results_[GROUP].tsv"
                ),
            )] = None,
    z_cutoff: Annotated[
        Optional[int],
        typer.Option(
            "-z",
            "--z_cutoff",
            "--zscore_cutoff",
            "--remove_outliers",
            "--remove_extreme",
            "--z_boundary",
            "--zcore_boundary",
            help=("Remove SNPs and associated regions whose counts exceed Z-Score cutoff.\n"
                  "Removing extreme outliers can provide extra layer of QC when measuring allelic imbalance. "
                  "(Default: None)"
                  )
            )
        ] = None,
):
    
    if len(groups) > 0:
        groups=groups[0]
    else:
        groups=None
    
    # Run single cell analysis
    run_ai_analysis_sc(count_file=counts,
                       bc_map=bc_map,
                       min_count=min,
                       pseudocount=pseudocount,
                       phase=phased,
                       sample=sample,
                       groups=groups,
                       out_file=out_file,
                       z_cutoff=z_cutoff
                       )


@app.command()
def compare_imbalance(
    counts: Annotated[str, typer.Argument(help="Count File")],
    bc_map: Annotated[str, typer.Argument(
        help=(
            "Two Column TSV file mapping specific barcodes to some grouping/celltype. "
            "Each line following format [BARCODE]\t[GROUP]"
            )
        )
                      ],
    min: Annotated[
        Optional[int],
        typer.Option(
            "--min",
            "--min_count",
            help=("Minimum allele count for measuring imbalance."
                  " (Default: 10)"
                  )
            )
        ] = None,
    pseudocount: Annotated[
        Optional[int],
        typer.Option(
            "-p",
            "--ps",
            "--pseudo",
            "--pseudocount",
            help=("Pseudocount added when measuring allelic imbalance. "
                  "(Default: 1)"
                  )
            )
        ] = None,
    sample: Annotated[
        Optional[str],
        typer.Option(
            "--sample",
            "--samp",
            "-s",
            help=(
                "Use heterozygous genotypes for this sample in count file. "
                "Automatically parses if data contains 0 or 1 sample. "
                "REQUIRED IF COUNT DATA CONTAINS MULTIPLE SAMPLES."
                ),
            )] = None,
    groups: Annotated[
        Optional[List[str]],
        typer.Option(
            "--groups",
            "--group",
            "--celltypes",
            "--g",
            help=(
                "Specific groups in barcode file/bc_map to compare allelic imbalance between. "
                "If providing input, requires a minimum of 2 groups. "
                "Uses all group combinations in barcode file/bc_map by default."
                ),
            )] = None,
    phased: Annotated[
        Optional[bool],
        typer.Option(
            "--phased/--unphased",
            help=(
                "If genotypes are phased use phasing information to measure imbalance.\n"
                "Otherwise or if --unphased selected, assume all haplotypes are equally likely during analysis.\n"
                "Autoparses genotype data by default if not denoted."
                ))] = None,
    out_file: Annotated[
        Optional[str],
        typer.Option(
            "--out_file",
            "--outfile",
            "--output",
            "--out",
            "-o",
            help=(
                "Output file for comparisons. "
                "Defaults to ai_results_[GROUP1]_[GROUP2].tsv"
                ),
            )] = None,
    z_cutoff: Annotated[
        Optional[int],
        typer.Option(
            "-z",
            "--z_cutoff",
            "--zscore_cutoff",
            "--remove_outliers",
            "--remove_extreme",
            "--z_boundary",
            "--zcore_boundary",
            help=("Remove SNPs and associated regions whose counts exceed Z-Score cutoff.\n"
                  "Removing extreme outliers can provide extra layer of QC when measuring allelic imbalance. "
                  "(Default: None)"
                  )
            )
        ] = None,
):
    
    if len(groups) > 0:
        groups=groups[0]
    else:
        groups=None
    
    # Run comparison
    run_ai_comparison(count_file=counts,
                      bc_map=bc_map,
                      min_count=min,
                      pseudocount=pseudocount,
                      phase=phased,
                      sample=sample,
                      groups=groups,
                      out_file=out_file,
                      z_cutoff=z_cutoff
                      )
    

if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()