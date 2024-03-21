from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from run_analysis import run_ai_analysis

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
    model: Annotated[
        Optional[str],
        typer.Option(
            "-m",
            "--model",
            help=(
                "Model used for measuring optimization parameter when finding imbalance. "
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
                    out_file=out_file,
                    region_col=region_col,
                    groupby=groupby)
    
    # TODO TEST CASES FOR TYPER


if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()