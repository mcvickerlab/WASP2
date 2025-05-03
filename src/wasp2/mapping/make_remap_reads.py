import timeit
import shutil
import tempfile
import subprocess
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple

import polars as pl
import numpy as np

import pysam
from pysam.libcalignmentfile import AlignmentFile

# Local imports
from wasp2.mapping.intersect_variant_data import make_intersect_df
from wasp2.mapping.remap_utils import (
    paired_read_gen, 
    paired_read_gen_stat, 
    get_read_het_data, 
    make_phased_seqs, 
    make_multi_seqs, 
    write_read
)


class ReadStats:
    """Track information about reads and SNPs that they overlap."""

    def __init__(self) -> None:
        # number of reads discarded because not mapped
        self.discard_unmapped: int = 0
        
        # number of reads discarded because not proper pair
        self.discard_improper_pair: int = 0

        # number of reads discarded because secondary match
        self.discard_secondary: int = 0

        # number of chimeric reads discarded
        self.discard_supplementary: int = 0

        # number of reads discarded because of too many allelic combinations
        self.discard_excess_reads: int = 0

        # number of reads discarded because expected pair is missing
        self.discard_missing_pair: int = 0
        
        # number of read pairs to remap
        self.remap_pair: int = 0
        
        # number of new pairs written
        self.write_pair: int = 0


def write_remap_bam(bam_file: str, intersect_file: str, r1_out: str, r2_out: str,
                    samples: List[str], max_seqs: int = 64) -> None:
    """
    Write remapped reads to FASTQ files after filtering and processing.

    This function first creates an intersection DataFrame using :func:`make_intersect_df` and then
    processes the reads from the BAM file according to variant intersection data. Depending on whether
    multiple samples are provided, it calls either :func:`swap_chrom_alleles_multi` or :func:`swap_chrom_alleles`
    to process each chromosome. Finally, it collates FASTQ files (for read1 and read2) from a temporary directory.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    intersect_file : str
        Path to the intersection file (from variant intersections).
    r1_out : str
        Path to the output FASTQ file for read1.
    r2_out : str
        Path to the output FASTQ file for read2.
    samples : list of str
        List of sample names.
    max_seqs : int, optional
        Maximum number of sequences to process (default is 64).

    Returns
    -------
    None

    Notes
    -----
    The following commented-out code is preserved:
    
        # TRY USING A CLASS OBJ
    """
    intersect_df = make_intersect_df(intersect_file, samples)
    
    # TRY USING A CLASS OBJ
    read_stats = ReadStats()
    
    # Use a temporary directory for intermediate chromosome FASTQ files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Previously, one could filter based on BAM header references:
        # remap_chroms = [c for c in bam.header.references if c in intersect_df.get_column("chrom").unique()]
        # Might need to change this/keep unordered for multiprocessed version
        remap_chroms = [c for c in intersect_df.get_column("chrom").unique(maintain_order=True)]

        if len(samples) > 1:
            for chrom in remap_chroms:
                swap_chrom_alleles_multi(
                    bam_file=bam_file, out_dir=tmpdir, df=intersect_df,
                    chrom=chrom, read_stats=read_stats
                )
        else:
            # Change from loop to multiprocess later
            for chrom in remap_chroms:
                swap_chrom_alleles(
                    bam_file=bam_file, out_dir=tmpdir, df=intersect_df,
                    chrom=chrom, read_stats=read_stats
                )
        
        # Get read1 FASTQ files from temporary directory
        r1_files = list(Path(tmpdir).glob("*_r1.fq"))
        with open(r1_out, "wb") as outfile_r1:
            for f in r1_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile_r1)
        
        # Get read2 FASTQ files
        r2_files = list(Path(tmpdir).glob("*_r2.fq"))
        with open(r2_out, "wb") as outfile_r2:
            for f in r2_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile_r2)
    
    print(f"Reads to remapped written to \n{r1_out}\n{r2_out}")


def swap_chrom_alleles(bam_file: str, out_dir: str, df: pl.DataFrame,
                       chrom: str, read_stats: ReadStats) -> None:
    """
    Swap alleles for a specific chromosome and write remapped reads.

    This function processes reads for a given chromosome by swapping alleles based on heterozygous
    SNP data. It partitions reads into dictionaries by their mate information, retrieves expected haplotype
    sequences using functions like :func:`get_read_het_data`, and writes new read pairs if allele swaps occur.
    After processing, it collates and writes FASTQ files for the chromosome.

    Parameters
    ----------
    bam_file : str
        Path to the original BAM file.
    out_dir : str
        Directory where the output files for the chromosome will be written.
    df : polars.DataFrame
        DataFrame containing intersection information for the variants.
    chrom : str
        Chromosome name to process.
    read_stats : ReadStats
        Object for tracking statistics of the remapping process.

    Returns
    -------
    None

    Notes
    -----
    """
    # Get hap columns (typically the last two columns)
    hap_cols = list(df.columns[-2:])
    # Create Chromosome-specific DataFrame (sort by start)
    chrom_df = df.filter(pl.col("chrom") == chrom).sort("start")
    
    r1_het_dict = chrom_df.filter(pl.col("mate") == 1).partition_by("read", as_dict=True, maintain_order=True)
    r2_het_dict = chrom_df.filter(pl.col("mate") == 2).partition_by("read", as_dict=True, maintain_order=True)
    
    # Create output BAM filename for the chromosome
    out_bam = str(Path(out_dir) / f"swapped_alleles_{chrom}.bam")
    
    start_chrom = timeit.default_timer()
    
    with AlignmentFile(bam_file, "rb") as bam, AlignmentFile(out_bam, "wb", header=bam.header) as out_file:
        if chrom not in bam.header.references:
            print(f"Skipping missing chrom: {chrom}")
            return
        
        for read1, read2 in paired_read_gen_stat(bam, read_stats, chrom=chrom):
            read_stats.remap_pair += 1
            og_name = read1.query_name
            r1_og_seq = read1.query_sequence
            r1_align_pos = read1.reference_start
            r2_og_seq = read2.query_sequence
            r2_align_pos = read2.reference_start
            
            write_num = 0  # Counter for reads written
            
            # Get SNP DataFrame for this read name
            r1_df = r1_het_dict.get(og_name)
            r2_df = r2_het_dict.get(og_name)
            
            # Og version using a function
            if r1_df is not None:
                r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
                if r1_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped position; skip read pair
                    continue
                r1_hap_list = [*make_phased_seqs(r1_het_data[0], *r1_het_data[1])]
            else:
                r1_hap_list = [r1_og_seq, r1_og_seq]
            
            if r2_df is not None:
                r2_het_data = get_read_het_data(r2_df, read2, hap_cols)
                if r2_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped position; skip read pair
                    continue
                r2_hap_list = [*make_phased_seqs(r2_het_data[0], *r2_het_data[1])]
            else:
                r2_hap_list = [r2_og_seq, r2_og_seq]
            
            # Create pairs to write only if at least one allele has changed
            write_pair_list = [
                (r1_hap_seq, r2_hap_seq)
                for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list)
                if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq)
            ]
            
            write_total = len(write_pair_list)
            
            # Write out new read pairs
            for r1_hap_seq, r2_hap_seq in write_pair_list:
                write_num += 1
                new_read_name = f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
                write_read(out_file, read1, r1_hap_seq, new_read_name)
                write_read(out_file, read2, r2_hap_seq, new_read_name)
                read_stats.write_pair += 1
            
        print(f"{chrom}: Processed in {timeit.default_timer() - start_chrom:.2f} seconds")
    
    # Collate and write out FASTQ files from the chromosome BAM file
    r1_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r1.fq")
    r2_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r2.fq")
    
    collate_cmd = ["samtools", "collate", "-u", "-O", out_bam]
    fastq_cmd = ["samtools", "fastq", "-1", r1_out, "-2", r2_out]
    
    collate_process = subprocess.run(collate_cmd, stdout=subprocess.PIPE, check=True)
    fastq_process = subprocess.run(fastq_cmd, input=collate_process.stdout, check=True)


def swap_chrom_alleles_multi(bam_file: str, out_dir: str, df: pl.DataFrame,
                             chrom: str, read_stats: ReadStats) -> None:
    """
    Swap alleles for a specific chromosome (multiple samples version) and write remapped reads.

    This function processes reads for a given chromosome by swapping alleles based on heterozygous
    SNP data across multiple samples. It partitions reads into dictionaries by their mate information, 
    extracts unique haplotype sequences using a combination of functions, and writes new read pairs if allele 
    swaps occur. Finally, it collates and writes FASTQ files for the processed chromosome.

    Parameters
    ----------
    bam_file : str
        Path to the original BAM file.
    out_dir : str
        Directory where the output files for the chromosome will be written.
    df : polars.DataFrame
        DataFrame containing intersection information for variants.
    chrom : str
        Chromosome name to process.
    read_stats : ReadStats
        Object for tracking statistics of the remapping process.

    Returns
    -------
    None
    """
    # Column data
    df_cols = df.columns[:5]
    hap_cols = df.columns[5:]
    
    # Create chromosome-specific DataFrame sorted by start position
    chrom_df = df.filter(pl.col("chrom") == chrom).sort("start")
    
    r1_het_dict = chrom_df.filter(pl.col("mate") == 1).partition_by("read", as_dict=True, maintain_order=True)
    r2_het_dict = chrom_df.filter(pl.col("mate") == 2).partition_by("read", as_dict=True, maintain_order=True)
    
    # Create output BAM filename for the chromosome
    out_bam = str(Path(out_dir) / f"swapped_alleles_{chrom}.bam")
    
    start_chrom = timeit.default_timer()
    
    with AlignmentFile(bam_file, "rb") as bam, AlignmentFile(out_bam, "wb", header=bam.header) as out_file:
        if chrom not in bam.header.references:
            print(f"Skipping missing chrom: {chrom}")
            return
        
        for read1, read2 in paired_read_gen_stat(bam, read_stats, chrom=chrom):
            read_stats.remap_pair += 1
            og_name = read1.query_name
            r1_og_seq = read1.query_sequence
            r1_align_pos = read1.reference_start
            r2_og_seq = read2.query_sequence
            r2_align_pos = read2.reference_start
            
            write_num = 0  # Counter that tracks reads written
            
            # Get SNP DataFrame for this read name from both dictionaries
            r1_df = r1_het_dict.get(og_name)
            r2_df = r2_het_dict.get(og_name)
            
            # Og version using a function
            if r1_df is not None:
                r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
                if r1_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped position; skip read pair
                    continue
                r1_hap_list = [*make_phased_seqs(r1_het_data[0], *r1_het_data[1])]
            else:
                r1_hap_list = [r1_og_seq, r1_og_seq]
            
            if r2_df is not None:
                r2_het_data = get_read_het_data(r2_df, read2, hap_cols)
                if r2_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped position; skip read pair
                    continue
                r2_hap_list = [*make_phased_seqs(r2_het_data[0], *r2_het_data[1])]
            else:
                r2_hap_list = [r2_og_seq, r2_og_seq]
            
            # Create pairs to write only if at least one allele is different from original
            write_pair_list = [
                (r1_hap_seq, r2_hap_seq)
                for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list)
                if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq)
            ]
            
            write_total = len(write_pair_list)
            
            # Write new read pairs to output BAM
            for r1_hap_seq, r2_hap_seq in write_pair_list:
                write_num += 1
                new_read_name = f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
                write_read(out_file, read1, r1_hap_seq, new_read_name)
                write_read(out_file, read2, r2_hap_seq, new_read_name)
                read_stats.write_pair += 1
            
        print(f"{chrom}: Processed in {timeit.default_timer() - start_chrom:.2f} seconds")
    
    # Collate and write out FASTQ files for the chromosome
    r1_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r1.fq")
    r2_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r2.fq")
    
    collate_cmd = ["samtools", "collate", "-u", "-O", out_bam]
    fastq_cmd = ["samtools", "fastq", "-1", r1_out, "-2", r2_out]
    
    collate_process = subprocess.run(collate_cmd, stdout=subprocess.PIPE, check=True)
    fastq_process = subprocess.run(fastq_cmd, input=collate_process.stdout, check=True)


# def swap_chrom_alleles(bam_file, out_dir, df, chrom, read_stats):
#     
#     # Get hap columns
#     hap_cols = list(df.columns[-2:])
#     # hap1_col, hap2_col = df.columns[-2:]
#     
#     # Create Chrom DF
#     
#     # Why is og order not maintained? Figure out and could skip sort step
#     chrom_df = df.filter(pl.col("chrom") == chrom).sort("start")
#     
#     r1_het_dict = chrom_df.filter(pl.col("mate") == 1).partition_by(
#         "read", as_dict=True, maintain_order=True)
#     
#     r2_het_dict = chrom_df.filter(pl.col("mate") == 2).partition_by(
#         "read", as_dict=True, maintain_order=True)
#     
#     # create chrom file
#     out_bam = str(Path(out_dir) / f"swapped_alleles_{chrom}.bam")
#     
#     # Might use to write per chrom stats later
#     # chrom_read_count = 0 
#     # chrom_write_count = 0
# 
#     start_chrom = timeit.default_timer()
#     
#     with AlignmentFile(bam_file, "rb") as bam, AlignmentFile(out_bam, "wb", header=bam.header) as out_file:
#         
#         if chrom not in bam.header.references:
#             print(f"Skipping missing chrom: {chrom}")
#             return
#         
#         for read1, read2 in paired_read_gen_stat(bam, read_stats, chrom=chrom):
#             
#             # chrom_read_count += 1
#             read_stats.remap_pair += 1
#             og_name = read1.query_name
#             r1_og_seq = read1.query_sequence
#             r1_align_pos = read1.reference_start
#             r2_og_seq = read2.query_sequence
#             r2_align_pos = read2.reference_start
#             
#             write_num = 0 # Counter that tracks reads written
#             
#             # Get snp df
#             r1_df = r1_het_dict.get(og_name)
#             r2_df = r2_het_dict.get(og_name)
#             
#             
#             # Og version using a func
#             if r1_df is not None:
#                 r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
#                 
#                 if r1_het_data is None:
#                     read_stats.discard_unmapped += 1
#                     # SNP overlaps unmapped pos
#                     continue
#                 r1_hap_list = [*make_phased_seqs(r1_het_data[0], *r1_het_data[1])]
#             
#             else:
#                 r1_hap_list = [r1_og_seq, r1_og_seq]
# 
#             
#             if r2_df is not None:
#                 r2_het_data = get_read_het_data(r2_df, read2, hap_cols)
#                 
#                 if r2_het_data is None:
#                     read_stats.discard_unmapped += 1
#                     # SNP overlaps unmapped pos
#                     continue
#                 
#                 r2_hap_list = [*make_phased_seqs(r2_het_data[0], *r2_het_data[1])]
# 
#             else:
#                 r2_hap_list = [r2_og_seq, r2_og_seq]
#             
#             # Create pairs to write
#             write_pair_list = [(r1_hap_seq, r2_hap_seq)
#                                for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list)
#                                if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq)]
#             
#             write_total = len(write_pair_list)
#             
#             # Get read pairs
#             for r1_hap_seq, r2_hap_seq in write_pair_list:
#                 write_num += 1
#                 new_read_name = f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
#                 write_read(out_file, read1, r1_hap_seq, new_read_name)
#                 write_read(out_file, read2, r2_hap_seq, new_read_name)
#                 read_stats.write_pair += 1
#                 # chrom_write_count += 1
# 
#         # WOWOW
#         # print(f"{chrom}: Processed {read_stats.remap_pair} pairs and wrote {read_stats.write_pair} new pairs in {timeit.default_timer() - start_chrom:.2f} seconds")
#     
#     # Collate and write out fastq now
#     collate_bam = str(Path(out_dir) / f"collate_{chrom}.bam")
#     r1_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r1.fq")
#     r2_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r2.fq")
#     
#     # Do I need to make another file???
#     pysam.collate(out_bam, "-o", collate_bam, catch_stdout=False)
#     pysam.fastq(collate_bam, "-1", r1_out, "-2", r2_out, catch_stdout=False)
#     # print(f"Created fastqs to be remapped in {Path(out_dir) / 'swapped_alleles_{chrom}_r*.fq'}")
