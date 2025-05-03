import timeit
import functools
import tempfile
import json
import warnings
from pathlib import Path
from typing import Any, Callable, Optional

# Import from local scripts
from wasp2.mapping.wasp_data_files import WaspDataFiles
from wasp2.mapping.intersect_variant_data import vcf_to_bed, process_bam, intersect_reads
from wasp2.mapping.make_remap_reads import write_remap_bam
from wasp2.mapping.filter_remap_reads import filt_remapped_reads, merge_filt_bam


# Decorator and Parser for read generation step
def tempdir_decorator(func: Callable) -> Callable:
    """Checks and creates a temporary directory for run_make_remap_reads().

    This decorator ensures that if the keyword argument 'temp_loc' is not provided when calling
    the decorated function, a temporary directory is created and passed as 'temp_loc'.

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
    def tempdir_wrapper(*args: Any, **kwargs: Any) -> Any:
        if kwargs.get("temp_loc", None) is not None:
            return func(*args, **kwargs)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                kwargs["temp_loc"] = tmpdir
                return func(*args, **kwargs)
    return tempdir_wrapper


@tempdir_decorator
def run_make_remap_reads(
    bam_file: str,
    vcf_file: str,
    is_paired: Optional[bool] = None,
    samples: Optional[Any] = None,
    is_phased: Optional[bool] = None,
    out_dir: Optional[str] = None,
    temp_loc: Optional[str] = None,
    out_json: Optional[str] = None
) -> None:
    """
    Parse initial inputs, find intersecting variants, and generate swapped allele reads for remapping.

    This function creates a WaspDataFiles object using the provided BAM and VCF files along with
    additional parameters. It then creates intermediary files by converting the VCF to BED, processing
    the BAM to extract remapped reads, and intersecting these with the variant BED file. Finally, it
    writes out FASTQ files containing the swapped allele reads and exports the file data as JSON.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    is_paired : optional bool
        Flag indicating whether the input BAM file is paired-end. Defaults to None.
    samples : optional
        Sample(s) to use for filtering. Defaults to None.
    is_phased : optional bool
        Flag indicating whether the VCF is phased. Defaults to None.
    out_dir : optional str
        Output directory for intermediary and final files. Defaults to None.
    temp_loc : optional str
        Directory for temporary files. Defaults to None.
    out_json : optional str
        Path to the output JSON file that will contain WASP file info. Defaults to None.

    Returns
    -------
    None
    """
    # Create Data Files
    wasp_files = WaspDataFiles(
        bam_file, vcf_file,
        is_paired=is_paired,
        samples=samples,
        is_phased=is_phased,
        out_dir=out_dir,
        temp_loc=temp_loc
    )
    
    # print(*vars(wasp_files).items(), sep="\n")
    
    # Create Checks for not integrated options
    if not wasp_files.is_paired:
        raise ValueError("Single-End not Implemented")
    
    if not wasp_files.is_phased:
        raise ValueError("Unphased not Implemented")
    
    if wasp_files.samples is None:
        raise ValueError("Zero samples not supported yet")
    
    # Should I create cache that checks for premade files??
    Path(wasp_files.out_dir).mkdir(parents=True, exist_ok=True)
    
    # Create Intermediary Files
    vcf_to_bed(
        vcf_file=wasp_files.vcf_file,
        out_bed=wasp_files.vcf_bed,
        samples=wasp_files.samples
    )

    process_bam(
        bam_file=wasp_files.bam_file,
        vcf_bed=wasp_files.vcf_bed,
        remap_bam=wasp_files.to_remap_bam,
        remap_reads=wasp_files.remap_reads,
        keep_bam=wasp_files.keep_bam,
        is_paired=wasp_files.is_paired
    )

    intersect_reads(
        remap_bam=wasp_files.to_remap_bam,
        vcf_bed=wasp_files.vcf_bed,
        out_bed=wasp_files.intersect_file
    )
    
    # print("INTERSECTION COMPLETE")
    
    # If a tempdir already exists??
    # Create remap fq
    write_remap_bam(
        wasp_files.to_remap_bam,
        wasp_files.intersect_file,
        wasp_files.remap_fq1,
        wasp_files.remap_fq2,
        wasp_files.samples
    )
    
    # print("WROTE READS TO BE REMAPPED")
    
    wasp_files.write_data(out_file=out_json)  # export JSON
    # print(f"File Data written to JSON...\n{out_json}")


# Decorator and Parser for post remap filtering
def check_filt_input(func: Callable) -> Callable:
    """Decorator that validates input types for run_wasp_filt().

    This decorator checks that either a JSON file containing WASP data is provided or that both
    'to_remap_bam' and 'keep_bam' are provided. It also sets a default output name if 'wasp_out_bam'
    is not given.

    Parameters
    ----------
    func : callable
        The function to decorate.

    Raises
    ------
    ValueError
        If neither wasp_data_json nor both to_remap_bam and keep_bam are provided.

    Returns
    -------
    callable
        The wrapped function with validated input.
    """
    @functools.wraps(func)
    def filt_wrapper(*args: Any, **kwargs: Any) -> Any:
        # Check if to_remap and keep bam given
        bam_input = all(
            (kwargs.get("to_remap_bam", None),
             kwargs.get("keep_bam", None))
        )
        
        # If JSON used instead of BAMs
        if kwargs.get("wasp_data_json", None) is not None:
            with open(kwargs.pop("wasp_data_json"), "r") as json_file:
                json_dict = json.load(json_file)
            
            if bam_input or any((kwargs.get("to_remap_bam", None),
                                 kwargs.get("keep_bam", None))):
                # Raise warning if both JSON and BAMs are provided
                warnings.warn(
                    ("Provided to_remap_bam+keep_bam ignored, using json input\n"
                     "Recommended Input of EITHER:\n"
                     "1. wasp_data_json\n2. to_remap_bam AND keep_bam")
                )
            
            # Set JSON inputs to BAMs
            kwargs["to_remap_bam"] = json_dict["to_remap_bam"]
            kwargs["keep_bam"] = json_dict["keep_bam"]
            
        elif not bam_input:
            raise ValueError(
                "Must provide either wasp_data_json OR BOTH to_remap_bam + keep_bam")
        
        elif "wasp_data_json" in kwargs:
            # Remove key if exists with None
            kwargs.pop("wasp_data_json")
        
        # Create default name if wasp_out_bam not given
        if kwargs.get("wasp_out_bam", None) is None:
            try:
                out_dir = json_dict["out_dir"]
                bam_prefix = json_dict["bam_prefix"]
            except Exception:
                out_dir = Path(kwargs["keep_bam"]).parent
                bam_prefix = Path(kwargs["keep_bam"]).name.rsplit("_keep.bam")[0]
            
            kwargs["wasp_out_bam"] = f"{out_dir}/{bam_prefix}_wasp_filt.bam"
            
        return func(*args, **kwargs)

    return filt_wrapper


@check_filt_input
def run_wasp_filt(
    remapped_bam: str,
    to_remap_bam: str,
    keep_bam: str,
    wasp_out_bam: str,
    remap_keep_bam: Optional[str] = None,
    remap_keep_file: Optional[str] = None
) -> None:
    """
    Filter remapped reads and merge with non-remapped reads to create a complete WASP filtered BAM.

    This function uses the provided remapped BAM and original BAM (to_remap_bam) along with a BAM
    containing reads that were not remapped (keep_bam) to filter out reads that did not remap to their
    expected location. It then merges the kept reads with the filtered remapped reads and outputs the final
    WASP filtered BAM file. Intermediate output can be optionally written to remap_keep_bam and remap_keep_file.

    Parameters
    ----------
    remapped_bam : str
        Path to the BAM file containing remapped reads.
    to_remap_bam : str
        Path to the original BAM file used to generate swapped alleles.
    keep_bam : str
        Path to the BAM file containing reads that were not remapped.
    wasp_out_bam : str
        Path to the output BAM file that will contain the final WASP filtered reads.
    remap_keep_bam : str, optional
        Path to an output BAM file for remapped reads that are kept. Defaults to None.
    remap_keep_file : str, optional
        Path to a TXT file with kept read names. Defaults to None.

    Returns
    -------
    None
    """
    # Handle temp
    if remap_keep_bam is None:
        with tempfile.TemporaryDirectory() as tmpdir:
            remap_keep_bam = f"{tmpdir}/wasp_remap_filt.bam"
            filt_remapped_reads(to_remap_bam, remapped_bam,
                                remap_keep_bam, keep_read_file=remap_keep_file)
            merge_filt_bam(keep_bam, remap_keep_bam, wasp_out_bam)
    else:
        filt_remapped_reads(to_remap_bam, remapped_bam, remap_keep_bam,
                            keep_read_file=remap_keep_file)
        print(f"\nWrote remapped bam with filtered reads to...\n{remap_keep_bam}\n")
        merge_filt_bam(keep_bam, remap_keep_bam, wasp_out_bam)
    
    # Finished
    print(f"\nWASP filtered Bam written to...\n{wasp_out_bam}\n")
