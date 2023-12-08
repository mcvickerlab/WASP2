from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from run_counting import run_count_variants

# app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
app = typer.Typer(pretty_exceptions_short=False)

# TODO GOTTA TEST THIS

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
                     )] = False
):
    
    # Parse sample string
    # print(samples)
    if len(samples) > 0:
        samples=samples[0]
    else:
        samples=None
    
    # print(samples)
    
    # run
    run_count_variants(bam_file=bam,
                       vcf_file=vcf,
                       region_file=region_file,
                       samples=samples,
                       use_region_names=use_region_names,
                       out_file=out_file,
                       temp_loc=temp_loc)
    
    # TODO TEST CASES FOR TYPER
    # TODO UNIT TEST NEXT



if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()