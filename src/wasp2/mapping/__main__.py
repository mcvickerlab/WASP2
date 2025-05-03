from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys

# Local Imports
from wasp2.mapping.run_mapping import run_make_remap_reads, run_wasp_filt

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
                "One or more samples to use in VCF "
                "Accepts comma delimited string, "
                "or file with one sample per line"
            )
        )
    ] = None,
    out_dir: Annotated[
        Optional[str],
        typer.Option(
            "--out_dir",
            "--outdir",
            "--out",
            "-o",
            help="Output directory for data to be remapped"
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
    out_json: Annotated[
        Optional[str],
        typer.Option(
            "--out_json",
            "--json",
            "--outjson",
            "-j",
            help=(
                "Output json containing WASP file info to this file instead of default. "
                "Defaults to [BAM_PREFIX]_wasp_data_files.json"
            )
        )
    ] = None,
    is_paired: Annotated[
        Optional[bool],
        typer.Option(
            "--paired/--single",
            help=(
                "Reads are paired or single. "
                "Will autoparse by default "
                "(SINGLE END NOT SUPPORTED YET)"
            )
        )
    ] = None,
    is_phased: Annotated[
        Optional[bool],
        typer.Option(
            "--phased/--unphased",
            help=(
                "If VCF is phased/unphased. "
                "Will autoparse by default "
                "(PHASED STRONGLY RECOMMENDED - SINGLE END NOT SUPPORTED YET)"
            )
        )
    ] = None,
) -> None:
    """
    Run the remapping pipeline to create reads for remapping.

    This command initiates the remapping pipeline by calling the 
    :func:`run_make_remap_reads` function with the provided BAM and VCF files, 
    along with optional sample filtering, output directory, temporary file location, 
    and JSON output file for WASP data.

    Parameters
    ----------
    bam : str
        Path to the BAM file.
    vcf : str
        Path to the VCF file.
    samples : Optional[List[str]]
        One or more samples to use in the VCF. Accepts comma delimited string or a file with one sample per line.
    out_dir : Optional[str]
        Output directory for data to be remapped.
    temp_loc : Optional[str]
        Directory for keeping intermediary files.
    out_json : Optional[str]
        Output JSON file for WASP file info. Defaults to "[BAM_PREFIX]_wasp_data_files.json" if not provided.
    is_paired : Optional[bool]
        Flag indicating whether the reads are paired. (SINGLE END NOT SUPPORTED YET)
    is_phased : Optional[bool]
        Flag indicating whether the VCF is phased or unphased.

    Returns
    -------
    None

    Examples
    --------
    >>> make_reads("input.bam", "variants.vcf", samples=["Sample1"], out_dir="output")
    
    Notes
    -----
    The following commented-out code is preserved:
    
        # print(samples)
        # print(samples)
    """
    # Parse sample string
    if samples is not None and len(samples) > 0:
        samples = samples[0]
    else:
        samples = None

    run_make_remap_reads(
        bam_file=bam,
        vcf_file=vcf,
        samples=samples,
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
        Optional[str],
        typer.Argument(
            help="to_remap_bam used to generate swapped alleles"
        )
    ] = None,
    keep_bam: Annotated[
        Optional[str],
        typer.Argument(
            help="BAM containing reads that were not remapped"
        )
    ] = None,
    wasp_data_json: Annotated[
        Optional[str],
        typer.Option(
            "--wasp_data_json",
            "--json",
            "-j",
            help=(
                "JSON containing WASP file info to load to_remap_bam and keep_bam"
            )
        )
    ] = None,
    out_bam: Annotated[
        Optional[str],
        typer.Option(
            "--out_bam",
            "--outbam",
            "--out",
            "-o",
            help=(
                "File to output filtered BAM. "
                "Will be created in default name and location if not provided"
            )
        )
    ] = None,
    remap_keep_bam: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_bam",
            help="Also output remapped BAM file containing kept reads"
        )
    ] = None,
    remap_keep_file: Annotated[
        Optional[str],
        typer.Option(
            "--remap_keep_file",
            help="Also output TXT file with kept read names"
        )
    ] = None
) -> None:
    """
    Filter the remapped BAM file to retain only reads that meet WASP filtering criteria.

    This command calls :func:`run_wasp_filt` with the provided remapped BAM file and optional parameters,
    including the original BAM file (to_remap_bam), a BAM file containing reads that were not remapped (keep_bam),
    and JSON output for WASP data.

    Parameters
    ----------
    remapped_bam : str
        Path to the remapped BAM file.
    to_remap_bam : Optional[str]
        Path to the original BAM file used to generate swapped alleles.
    keep_bam : Optional[str]
        BAM file containing reads that were not remapped.
    wasp_data_json : Optional[str]
        JSON file containing WASP file info to load to_remap_bam and keep_bam.
    out_bam : Optional[str]
        File to output the filtered BAM. If not provided, a default name and location are used.
    remap_keep_bam : Optional[str]
        Also output remapped BAM file containing kept reads.
    remap_keep_file : Optional[str]
        Also output a TXT file with kept read names.

    Returns
    -------
    None

    Examples
    --------
    >>> filter_remapped("remapped.bam", to_remap_bam="original.bam", keep_bam="not_remapped.bam", out_bam="filtered.bam")
    
    Notes
    -----
    The following commented-out code is preserved:
    
        # print(remapped_bam)
        # print(to_remap_bam)
        # print(keep_bam)
        # print(wasp_data_json)
        # print(out_bam)
        # print(remap_keep_bam)
        # print(remap_keep_file)
    """
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
