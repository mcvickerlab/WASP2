import tempfile
from pathlib import Path
import timeit
from typing import Optional

import pysam
from pysam.libcalignmentfile import AlignmentFile

from remap_utils import paired_read_gen

def filt_remapped_reads(
    to_remap_bam: str,
    remapped_bam: str,
    filt_out_bam: str,
    keep_read_file: Optional[str] = None
) -> None:
    
    pos_dict = {}
    total_dict = {}
    keep_set = set()
    
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
    
    start_time = timeit.default_timer()
    
    # Merge for for complete filt bam
    pysam.merge("-f", "-o", out_bam, keep_bam, remapped_filt_bam, catch_stdout=False)
    print(f"Merged BAM in {timeit.default_timer() - start_time:.2f} seconds")
    
    start_sort = timeit.default_timer()
    pysam.sort(out_bam, "-o", out_bam, catch_stdout=False)
    pysam.index(out_bam, catch_stdout=False)
    
    print(f"Sorted and Indexed BAM in {timeit.default_timer() - start_sort:.2f} seconds")
    
    # print(f"\nWrote merged WASP filtered BAM to...\n{out_bam}")