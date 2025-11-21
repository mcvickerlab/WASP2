import timeit
import functools
import tempfile
import json
import warnings
from pathlib import Path
from typing import Optional, Union, List, Callable, Any

# Import from local scripts
from .wasp_data_files import WaspDataFiles
from .intersect_variant_data import vcf_to_bed, process_bam, intersect_reads

from .make_remap_reads import write_remap_bam
from .filter_remap_reads import filt_remapped_reads, merge_filt_bam


# Decorator and Parser for read generation step
def tempdir_decorator(func: Callable[..., Any]) -> Callable[..., Any]:
    """Checks and makes tempdir for
    run_make_remap_reads()
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
    samples: Optional[Union[str, List[str]]] = None,
    is_phased: Optional[bool] = None,
    out_dir: Optional[str] = None,
    temp_loc: Optional[str] = None,
    out_json: Optional[str] = None
) -> None:
    """
    Parser that parses initial input.
    Finds intersecting variants and generates
    swapped allele reads to be remapped.
    

    :param bam_file: _description_
    :type bam_file: _type_
    :param vcf_file: _description_
    :type vcf_file: _type_
    :param is_paired: _description_, defaults to None
    :type is_paired: _type_, optional
    :param samples: _description_, defaults to None
    :type samples: _type_, optional
    :param is_phased: _description_, defaults to None
    :type is_phased: _type_, optional
    :param out_dir: _description_, defaults to None
    :type out_dir: _type_, optional
    :param temp_loc: _description_, defaults to None
    :type temp_loc: _type_, optional
    :param out_json: _description_, defaults to None
    :type out_json: _type_, optional
    """
    
    
    # Create Data Files
    wasp_files = WaspDataFiles(bam_file, vcf_file,
                               is_paired=is_paired,
                               samples=samples,
                               is_phased=is_phased,
                               out_dir=out_dir,
                               temp_loc=temp_loc)
    
    # print(*vars(wasp_files).items(), sep="\n")
    
    # Create Checks for not integrated options
    if not wasp_files.is_paired:
        raise ValueError("Single-End not Implemented")

    if not wasp_files.is_phased:
        raise ValueError("Unphased not Implemented")

    if wasp_files.samples is None:
        raise ValueError("Zero samples not supported yet")

    # Type narrowing: help mypy understand the types after the above checks
    # - is_paired is True, so remap_fq2 is str (not None)
    # - samples is List[str] (normalized in WaspDataFiles, not None)
    assert isinstance(wasp_files.samples, list), "samples should be normalized to list"
    assert wasp_files.remap_fq2 is not None, "remap_fq2 should be set when is_paired is True"

    # Should I create cache that checks for premade files??
    Path(str(wasp_files.out_dir)).mkdir(parents=True, exist_ok=True)


    # Create Intermediary Files
    vcf_to_bed(vcf_file=str(wasp_files.vcf_file),
               out_bed=wasp_files.vcf_bed,
               samples=wasp_files.samples)


    process_bam(bam_file=str(wasp_files.bam_file),
                vcf_bed=wasp_files.vcf_bed,
                remap_bam=wasp_files.to_remap_bam,
                remap_reads=wasp_files.remap_reads,
                keep_bam=wasp_files.keep_bam,
                is_paired=wasp_files.is_paired)


    intersect_reads(remap_bam=wasp_files.to_remap_bam,
                    vcf_bed=wasp_files.vcf_bed,
                    out_bed=wasp_files.intersect_file)
    
    
    # print("INTERSECTION COMPLETE")
    
    # If a tempdir already exists??

    # Create remap fq
    write_remap_bam(wasp_files.to_remap_bam,
                    wasp_files.intersect_file,
                    wasp_files.remap_fq1,
                    wasp_files.remap_fq2,
                    wasp_files.samples)
    
    
    # print("WROTE READS TO BE REMAPPED")
    
    
    wasp_files.write_data(out_file=out_json) # export json
    # print(f"File Data written to JSON...\n{out_json}")


# Decorator and Parser for post remap filtering
def check_filt_input(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator that parses valid input types
    for run_wasp_filt()

    :param func: _description_
    :type func: _type_
    :raises ValueError: _description_
    :return: _description_
    :rtype: _type_
    """

    @functools.wraps(func)
    def filt_wrapper(*args: Any, **kwargs: Any) -> Any:

        # Check if to_remap and keep bam given
        bam_input = all(
            (kwargs.get("to_remap_bam", None),
             kwargs.get("keep_bam", None))
        )
        
        # If json used instead of bams
        if kwargs.get("wasp_data_json", None) is not None:
            
            with open(kwargs.pop("wasp_data_json"), "r") as json_file:
                json_dict = json.load(json_file)
            
            if bam_input or any((kwargs.get("to_remap_bam", None),
                                 kwargs.get("keep_bam", None))):
                
                # Raise warning if json and bams given
                warnings.warn(
                    ("Provided to_remap_bam+keep_bam ignored, using json input\n"
                     "Recommended Input of EITHER:\n"
                     "1. wasp_data_json\n2. to_remap_bam AND keep_bam")
                )
            
            # Set json inputs to bams
            kwargs["to_remap_bam"] = json_dict["to_remap_bam"]
            kwargs["keep_bam"] = json_dict["keep_bam"]
            
        elif not bam_input:
            raise ValueError(
                "Must provide either wasp_data_json OR BOTH to_remap_bam + keep_bam")
        
        elif "wasp_data_json" in kwargs:
            # remove if None, but key exists in kwargs
            kwargs.pop("wasp_data_json")
            
        
        # Create default name if wasp_out_bam not given
        if kwargs.get("wasp_out_bam", None) is None:
            
            # If data included in json
            try:
                out_dir = json_dict["out_dir"]
                bam_prefix = json_dict["bam_prefix"]
            except:
                out_dir = Path(kwargs["keep_bam"]).parent
                bam_prefix = Path(kwargs["keep_bam"]).name.rsplit("_keep.bam")[0]
            
            # create output file
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
    Filter reads that remap to the same loc
    and merges with non-remapped reads to create
    a complete WASP filtered BAM

    :param remapped_bam: _description_
    :type remapped_bam: _type_
    :param to_remap_bam: _description_
    :type to_remap_bam: _type_
    :param keep_bam: _description_
    :type keep_bam: _type_
    :param wasp_out_bam: _description_
    :type wasp_out_bam: _type_
    :param remap_keep_bam: _description_, defaults to None
    :type remap_keep_bam: _type_, optional
    :param remap_keep_file: _description_, defaults to None
    :type remap_keep_file: _type_, optional
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
