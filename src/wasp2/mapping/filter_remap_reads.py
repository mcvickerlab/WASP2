import sys
import timeit
import re
import functools
import tempfile
import warnings
from pathlib import Path
from typing import Any, Dict, Set, Optional, Tuple

import pysam
from pysam.libcalignmentfile import AlignmentFile

from wasp2.mapping.remap_utils import paired_read_gen


def filt_remapped_reads(
    to_remap_bam: str,
    remapped_bam: str,
    filt_out_bam: str,
    keep_read_file: Optional[str] = None,
) -> None:
    """
    Filter remapped reads based on original mapping positions.

    This function processes a BAM file containing remapped reads by comparing the remapped
    positions with the expected original positions (encoded in the read names). It maintains
    dictionaries for the expected positions and total read counts and generates a set of read
    names that meet the criteria. Reads failing the criteria are removed. Finally, it filters
    the original BAM file to include only the reads in the "keep" set using pysam.view.

    Parameters
    ----------
    to_remap_bam : str
        Path to the original BAM file that was subject to remapping.
    remapped_bam : str
        Path to the BAM file containing remapped reads.
    filt_out_bam : str
        Path where the filtered BAM file will be written.
    keep_read_file : str, optional
        Path to a file where the read names to keep will be written. If None, a temporary file is used.

    Returns
    -------
    None

    Notes
    -----
    The following commented-out code is preserved:
    
        # print(f"{len(keep_set)} pairs remapped successfuly!")
        # print(f"{num_removed} pairs removed!") # Inaccurate?
        # print(vars(read_stats))
        # print(f"Wrote reads that successfully remapped to {keep_read_file}")
        # print(f"Wrote bam with filtered reads to {filt_out_bam}")
    """
    pos_dict: Dict[str, Tuple[int, int]] = {}
    total_dict: Dict[str, int] = {}
    keep_set: Set[str] = set()
    
    num_removed = 0
    
    with AlignmentFile(remapped_bam, "rb") as bam:
        # nostat???
        for read1, read2 in paired_read_gen(bam):
            read_name_split = read1.query_name.split("_WASP_")
            read_name = read_name_split[0]
            
            if read_name not in pos_dict:
                # First time seeing read, add to dict and set
                read_data = tuple(map(int, read_name_split[1].split("_", maxsplit=3)))
                pos_dict[read_name] = (read_data[0], read_data[1])
                total_dict[read_name] = read_data[3]
                keep_set.add(read_name)
            elif read_name not in keep_set:
                # If seen, but removed from set, skip
                # print(f"Removed {read_name} skipping {read1.query_name}")
                continue
            
            # Count down reads seen
            total_dict[read_name] -= 1
            
            # Check for equality
            if (read1.reference_start, read1.next_reference_start) != pos_dict[read_name]:
                keep_set.remove(read_name)
                total_dict.pop(read_name)
                num_removed += 1
            elif total_dict[read_name] == 0:
                # Found expected number of reads
                total_dict.pop(read_name)
                pos_dict.pop(read_name)
    
    # Remove reads with Missing Counts
    missing_count_set = set(total_dict.keys())
    num_removed += len(missing_count_set)
    keep_set = keep_set - missing_count_set

    # Write keep reads to file
    # print(f"{len(keep_set)} pairs remapped successfuly!")
    # print(f"{num_removed} pairs removed!") # Inaccurate?
    # print(vars(read_stats))
    
    # print(f"Wrote reads that successfully remapped to {keep_read_file}")
    
    # Check if need to create temp file
    if keep_read_file is None:
        with tempfile.NamedTemporaryFile("w") as file:
            file.write("\n".join(keep_set))
            pysam.view("-N", file.name, "-o", filt_out_bam, to_remap_bam, catch_stdout=False)
    else:
        with open(keep_read_file, "w") as file:
            file.write("\n".join(keep_set))
        
        print(f"\nWrote Remapped Reads kept to...\n{keep_read_file}\n")
        pysam.view("-N", keep_read_file, "-o", filt_out_bam, to_remap_bam, catch_stdout=False)
    
    # print(f"Wrote bam with filtered reads to {filt_out_bam}")


def merge_filt_bam(keep_bam: str, remapped_filt_bam: str, out_bam: str) -> None:
    """
    Merge and sort filtered BAM files.

    This function merges a BAM file containing kept reads and a BAM file with filtered remapped reads
    using pysam.merge, then sorts and indexes the resulting BAM file.

    Parameters
    ----------
    keep_bam : str
        Path to the BAM file containing reads to keep.
    remapped_filt_bam : str
        Path to the BAM file containing filtered remapped reads.
    out_bam : str
        Path where the merged, sorted, and indexed BAM file will be written.

    Returns
    -------
    None

    Notes
    -----
    The following commented-out code is preserved:
    
        # print(f"Wrote merged WASP filtered BAM to {out_bam}")
    """
    start_time = timeit.default_timer()
    
    # Merge for complete filt bam
    pysam.merge("-f", "-o", out_bam, keep_bam, remapped_filt_bam, catch_stdout=False)
    print(f"Merged BAM in {timeit.default_timer() - start_time:.2f} seconds")
    
    start_sort = timeit.default_timer()
    pysam.sort(out_bam, "-o", out_bam, catch_stdout=False)
    pysam.index(out_bam, catch_stdout=False)
    
    print(f"Sorted and Indexed BAM in {timeit.default_timer() - start_sort:.2f} seconds")
    
    # print(f"Wrote merged WASP filtered BAM to {out_bam}")
