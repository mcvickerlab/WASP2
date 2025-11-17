from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from .run_mapping import run_make_remap_reads, run_wasp_filt


app = typer.Typer()
# app = typer.Typer(pretty_exceptions_show_locals=False)
# app = typer.Typer(pretty_exceptions_short=False)

@app.command()
def make_reads(
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
                "One or more samples to use in VCF"
                "Accepts comma delimited string, "
                "or file with one sample per line"
            )
            )] = None,
    out_dir: Annotated[
        Optional[str],
        typer.Option(
            "--out_dir",
            "--outdir",
            "--out",
            "-o",
            help="Output directory for data to be remapped"
            )] = None,
    temp_loc: Annotated[
        Optional[str],
        typer.Option(
            "--temp_loc",
            "--temp",
            "-t",
            help=(
                "Directory for keeping intermediary files."
                "Defaults to removing intermediary files using temp directory")
            )] = None,
    # is_phased: Annotated[Optional[bool], typer.Argument()] = None,
    out_json: Annotated[
        Optional[str],
        typer.Option("--out_json",
                     "--json",
                     "--outjson",
                     "-j",
                     help=(
                         "Output json containing wasp file info to this file instead of default. "
                         "Defaults to [BAM_PREFIX]_wasp_data_files.json"
                     )
                     )] = None,
    is_paired: Annotated[
        Optional[bool],
        typer.Option("--paired/--single",
                     help=(
                         "Reads are paired or single. "
                         "Will autoparse by default "
                         "(SINGLE END NOT SUPPORTED YET)"
                         )
                     )] = None,
    is_phased: Annotated[
        Optional[bool],
        typer.Option("--phased/--unphased",
                     help=(
                         "If VCF is phased/unphased"
                         "Will autoparse by default "
                         "(PHASED STRONGLY RECOMMENDED-SINGLE END NOT SUPPORTED YET)"
                         )
                     )] = None,
) -> None:
    """Generate reads with swapped alleles for remapping."""

    # Parse sample string
    sample_str: Optional[str]
    if samples is not None and len(samples) > 0:
        sample_str = samples[0]
    else:
        sample_str = None

    run_make_remap_reads(
        bam_file=bam,
        vcf_file=vcf,
        samples=sample_str,
        out_dir=out_dir,
        temp_loc=temp_loc,
        out_json=out_json,
        is_paired=is_paired,
        is_phased=is_phased
        )




@app.command()
def filter_remapped(
    remapped_bam: Annotated[str, typer.Argument(help="remapped BAM File")],
    to_remap_bam: Annotated[
        Optional[str], typer.Argument(
            help="to_remap_bam used to generate swapped alleles")
        ] = None,
    keep_bam: Annotated[
        Optional[str], typer.Argument(
            help="BAM containing reads that were not remapped")
        ] = None,
    wasp_data_json: Annotated[
        Optional[str],
        typer.Option("--wasp_data_json",
                     "--json",
                     "-j",
                     help=(
                         "json containing wasp file info to load to_remap_bam and keep_bam"
                         ))
        ] = None,
    out_bam: Annotated[
        Optional[str],
        typer.Option(
            "--out_bam",
            "--outbam",
            "--out",
            "-o",
            help=(
                "File to output filt bam. "
                "Will be created in default name and loc if not provided"
                )
            )] = None,
    remap_keep_bam: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_bam",
            help=(
                "Also output remapped bam file containing kept reads"
                )
            )] = None,
    remap_keep_file: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_file",
            help=(
                "Also output txt file with kept read names"
                )
            )] = None
) -> None:
    """Filter remapped reads using WASP algorithm."""

    # Checks
    # print(remapped_bam)
    # print(to_remap_bam)
    # print(keep_bam)
    # print(wasp_data_json)
    # print(out_bam)
    # print(remap_keep_bam)
    # print(remap_keep_file)
    
    # Run WASP Filt
    run_wasp_filt(
        remapped_bam,
        to_remap_bam=to_remap_bam,
        keep_bam=keep_bam,
        wasp_out_bam=out_bam,
        remap_keep_bam=remap_keep_bam,
        remap_keep_file=remap_keep_file,
        wasp_data_json=wasp_data_json
        )
    

if __name__ == "__main__":
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()
