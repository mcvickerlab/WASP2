import timeit

import shutil
import tempfile

from pathlib import Path
from typing import List

# from collections import defaultdict

import polars as pl

import pysam
from pysam.libcalignmentfile import AlignmentFile

# local imports
from intersect_variant_data import make_intersect_df
from remap_utils import paired_read_gen, paired_read_gen_stat, get_read_het_data, make_phased_seqs, make_multi_seqs, write_read


# TRY subprocess
import subprocess


class ReadStats(object):
    """Track information about reads and SNPs that they overlap"""

    def __init__(self) -> None:
        # number of read matches to reference allele
        # self.ref_count = 0
        # number of read matches to alternative allele
        # self.alt_count = 0
        # number of reads that overlap SNP but match neither allele
        # self.other_count = 0

        # number of reads discarded becaused not mapped
        self.discard_unmapped = 0
        
        # number of reads discarded because not proper pair
        self.discard_improper_pair = 0

        # number of reads discarded because mate unmapped
        # self.discard_mate_unmapped = 0

        # paired reads map to different chromosomes
        # self.discard_different_chromosome = 0

        # number of reads discarded because secondary match
        self.discard_secondary = 0

        # number of chimeric reads discarded
        self.discard_supplementary = 0

        # number of reads discarded because of too many overlapping SNPs
        # self.discard_excess_snps = 0
        
        # number of reads discarded because too many allelic combinations
        self.discard_excess_reads = 0

        # when read pairs share SNP locations but have different alleles there
        # self.discard_discordant_shared_snp = 0
        
        # reads where we expected to see other pair, but it was missing
        # possibly due to read-pairs with different names
        self.discard_missing_pair = 0
        
        # number of single reads that need remapping
        # self.remap_single = 0
        
        # number of read pairs to remap
        self.remap_pair = 0
        
        # Number of new pairs written
        self.write_pair = 0


def write_remap_bam(
    bam_file: str,
    intersect_file: str,
    r1_out: str,
    r2_out: str,
    samples: List[str],
    max_seqs: int = 64
) -> None:
    intersect_df = make_intersect_df(intersect_file, samples)
    
    # TRY USING A CLASS OBJ
    read_stats = ReadStats()
    
    # Should this be r or rb? Need to figure out Errno 9 bad file descrip error
    # with AlignmentFile(bam_file, "rb") as bam, tempfile.TemporaryDirectory() as tmpdir:
    with tempfile.TemporaryDirectory() as tmpdir:
        
        # remap_chroms = [c for c in bam.header.references
        #                 if c in intersect_df.get_column("chrom").unique()]
        
        # Might need to change this/keep unordered for multiprocesed version
        remap_chroms = [c for c in intersect_df.get_column("chrom").unique(maintain_order=True)]

        if len(samples) > 1:
            for chrom in remap_chroms:
                swap_chrom_alleles_multi(bam_file=bam_file, out_dir=tmpdir,
                                         df=intersect_df, chrom=chrom,
                                         read_stats=read_stats)

        else:
            # tmpdir="/iblm/netapp/home/aho/projects/wasp/testing/mapping_v2/outputs/test_remap_v1/samp_cli_v1/chrom_files"

            # Change from loop to multiprocess later
            for chrom in remap_chroms:
                
                swap_chrom_alleles(bam_file=bam_file, out_dir=tmpdir,
                                   df=intersect_df, chrom=chrom,
                                   read_stats=read_stats)
        
        # Get r1 files
        r1_files = list(Path(tmpdir).glob("*_r1.fq"))
        
        with open(r1_out, "wb") as outfile_r1:
            for f in r1_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile_r1)
        
        
        r2_files = list(Path(tmpdir).glob("*_r2.fq"))
        
        with open(r2_out, "wb") as outfile_r2:
            for f in r2_files:
                with open(f, "rb") as infile:
                    shutil.copyfileobj(infile, outfile_r2)
    
    print(f"Reads to remapped written to \n{r1_out}\n{r2_out}")


def swap_chrom_alleles(
    bam_file: str,
    out_dir: str,
    df: pl.DataFrame,
    chrom: str,
    read_stats: ReadStats
) -> None:

    # Get hap columns
    hap_cols = list(df.columns[-2:])
    # hap1_col, hap2_col = df.columns[-2:]
    
    # Create Chrom DF
    
    # Why is og order not maintained? Figure out and could skip sort step
    chrom_df = df.filter(pl.col("chrom") == chrom).sort("start")

    # Polars 1.x: partition_by(as_dict=True) returns dict with tuple keys
    # Need to extract first element of tuple for single-column partition
    r1_partition = chrom_df.filter(pl.col("mate") == 1).partition_by(
        "read", as_dict=True, maintain_order=True)
    r1_het_dict = {k[0] if isinstance(k, tuple) else k: v for k, v in r1_partition.items()}

    r2_partition = chrom_df.filter(pl.col("mate") == 2).partition_by(
        "read", as_dict=True, maintain_order=True)
    r2_het_dict = {k[0] if isinstance(k, tuple) else k: v for k, v in r2_partition.items()}
    
    # create chrom file
    out_bam = str(Path(out_dir) / f"swapped_alleles_{chrom}.bam")
    
    # Might use to write per chrom stats later
    # chrom_read_count = 0 
    # chrom_write_count = 0

    start_chrom = timeit.default_timer()
    
    # Maybe check if file descrip not closed properly???
    with AlignmentFile(bam_file, "rb") as bam, AlignmentFile(out_bam, "wb", header=bam.header) as out_file:

        if chrom not in bam.header.references:
            print(f"Skipping missing chrom: {chrom}")
            return
        
        for read1, read2 in paired_read_gen_stat(bam, read_stats, chrom=chrom):
            
            # chrom_read_count += 1
            read_stats.remap_pair += 1
            og_name = read1.query_name
            r1_og_seq = read1.query_sequence
            r1_align_pos = read1.reference_start
            r2_og_seq = read2.query_sequence
            r2_align_pos = read2.reference_start
            
            write_num = 0 # Counter that tracks reads written
            
            # Get snp df
            r1_df = r1_het_dict.get(og_name)
            r2_df = r2_het_dict.get(og_name)
            
            
            # Og version using a func
            if r1_df is not None:
                r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
                
                if r1_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped pos
                    continue
                r1_hap_list = [*make_phased_seqs(r1_het_data[0], *r1_het_data[1])]
            
            else:
                r1_hap_list = [r1_og_seq, r1_og_seq]

            
            if r2_df is not None:
                r2_het_data = get_read_het_data(r2_df, read2, hap_cols)
                
                if r2_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped pos
                    continue
                
                r2_hap_list = [*make_phased_seqs(r2_het_data[0], *r2_het_data[1])]

            else:
                r2_hap_list = [r2_og_seq, r2_og_seq]
            
            # Create pairs to write
            write_pair_list = [(r1_hap_seq, r2_hap_seq)
                               for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list)
                               if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq)]
            
            write_total = len(write_pair_list)
            
            # Get read pairs
            for r1_hap_seq, r2_hap_seq in write_pair_list:
                write_num += 1
                new_read_name = f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
                write_read(out_file, read1, r1_hap_seq, new_read_name)
                write_read(out_file, read2, r2_hap_seq, new_read_name)
                read_stats.write_pair += 1
                # chrom_write_count += 1

        # print(f"{chrom}: Processed {read_stats.remap_pair} pairs and wrote {read_stats.write_pair} new pairs in {timeit.default_timer() - start_chrom:.2f} seconds")
        print(f"{chrom}: Processed in {timeit.default_timer() - start_chrom:.2f} seconds")

    # Collate and write out fastq
    r1_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r1.fq")
    r2_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r2.fq")
    
    # Do I need to make another file???
    
    # pysam.collate("-u","-o", collate_bam, out_bam, catch_stdout=False)
    # pysam.fastq("-1", r1_out, "-2", r2_out, collate_bam,
    #             "--verbosity", "0", catch_stdout=False)
    
    
    # TRY SUBPROCESS METHOD
    
    # TRY piping subprocess, so no pysam wrapper
    collate_cmd = ["samtools", "collate",
                   "-u", "-O", out_bam]
    
    fastq_cmd = ["samtools", "fastq",
                 "-1", r1_out, "-2", r2_out]
    
    collate_process = subprocess.run(collate_cmd, stdout=subprocess.PIPE, check=True)
    fastq_process = subprocess.run(fastq_cmd, input=collate_process.stdout, check=True)


def swap_chrom_alleles_multi(
    bam_file: str,
    out_dir: str,
    df: pl.DataFrame,
    chrom: str,
    read_stats: ReadStats
) -> None:

    # column data
    df_cols = df.columns[:5]
    hap_cols = df.columns[5:]
    
    # Create chrom df
    chrom_df = df.filter(pl.col("chrom") == chrom).sort("start")

    # Polars 1.x: partition_by(as_dict=True) returns dict with tuple keys
    r1_partition = chrom_df.filter(pl.col("mate") == 1).partition_by(
        "read", as_dict=True, maintain_order=True)
    r1_het_dict = {k[0] if isinstance(k, tuple) else k: v for k, v in r1_partition.items()}

    r2_partition = chrom_df.filter(pl.col("mate") == 2).partition_by(
        "read", as_dict=True, maintain_order=True)
    r2_het_dict = {k[0] if isinstance(k, tuple) else k: v for k, v in r2_partition.items()}
    
    
    # create chrom file
    out_bam = str(Path(out_dir) / f"swapped_alleles_{chrom}.bam") # temp, create correct in file data
    

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

            write_num = 0 # Counter that tracks reads written

            # Get snp_df
            r1_df = r1_het_dict.pop(og_name, None)
            r2_df = r2_het_dict.pop(og_name, None)
            
            if (r1_df is not None) and (r2_df is not None):
                read_df = r1_df.vstack(r2_df) # Combine for testing equality
            elif r1_df is not None:
                read_df = r1_df
            elif r2_df is not None:
                read_df = r2_df
            else:
                # TEMPORARY FIX FOR BUG????
                # NOT SURE WHY SOME READS WOULD SHOW UP BUT NOT OVERLAP A SNP
                continue


            # if (r1_df is not None) and (r2_df is not None):
            #     read_df = r1_df.vstack(r2_df) # Combine for testing equality
            # elif r1_df is not None:
            #     read_df = r1_df
            # else:
            #     read_df = r2_df


            # Get unique haps
            unique_cols = (
                read_df.select(
                    pl.col(hap_cols).str.concat("")
                ).transpose(
                    include_header=True, column_names=["hap"]
                ).unique(
                    subset=["hap"]).get_column("column")
            )


            # create new col data
            use_cols = [*df_cols, *unique_cols]
            num_haps = len(unique_cols)


            if r1_df is not None:
                r1_df = r1_df.select(pl.col(use_cols))

                r1_het_data = get_read_het_data(r1_df, read1, unique_cols)

                if r1_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped pos
                    continue

                r1_hap_list = make_multi_seqs(*r1_het_data)
            else:
                r1_hap_list = [r1_og_seq] * num_haps


            if r2_df is not None:
                r2_df = r2_df.select(pl.col(use_cols))

                r2_het_data = get_read_het_data(r2_df, read2, unique_cols)

                if r2_het_data is None:
                    read_stats.discard_unmapped += 1
                    # SNP overlaps unmapped pos
                    continue

                r2_hap_list = make_multi_seqs(*r2_het_data)
            else:
                r2_hap_list = [r2_og_seq] * num_haps



            # Create Pairs to write
            write_pair_list = [(r1_hap_seq, r2_hap_seq) 
                               for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list) 
                               if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq)]

            write_total = len(write_pair_list)

            # Get read pairs
            for r1_hap_seq, r2_hap_seq in write_pair_list:
                write_num += 1
                new_read_name = f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
                
                write_read(out_file, read1, r1_hap_seq, new_read_name)
                write_read(out_file, read2, r2_hap_seq, new_read_name)
                read_stats.write_pair += 1
        
        # Done
        print(f"{chrom}: Processed in {timeit.default_timer() - start_chrom:.2f} seconds")    

    # Collate and write out fastq
    r1_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r1.fq")
    r2_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r2.fq")
    
    collate_cmd = ["samtools", "collate",
                   "-u", "-O", out_bam]
    
    fastq_cmd = ["samtools", "fastq",
                 "-1", r1_out, "-2", r2_out]
    
    collate_process = subprocess.run(collate_cmd, stdout=subprocess.PIPE, check=True)
    fastq_process = subprocess.run(fastq_cmd, input=collate_process.stdout, check=True)





# def swap_chrom_alleles(bam_file, out_dir, df, chrom, read_stats):
    
#     # Get hap columns
#     hap_cols = list(df.columns[-2:])
#     # hap1_col, hap2_col = df.columns[-2:]
    
#     # Create Chrom DF
    
#     # Why is og order not maintained? Figure out and could skip sort step
#     chrom_df = df.filter(pl.col("chrom") == chrom).sort("start")
    
#     r1_het_dict = chrom_df.filter(pl.col("mate") == 1).partition_by(
#         "read", as_dict=True, maintain_order=True)
    
#     r2_het_dict = chrom_df.filter(pl.col("mate") == 2).partition_by(
#         "read", as_dict=True, maintain_order=True)
    
#     # create chrom file
#     out_bam = str(Path(out_dir) / f"swapped_alleles_{chrom}.bam")
    
#     # Might use to write per chrom stats later
#     # chrom_read_count = 0 
#     # chrom_write_count = 0

#     start_chrom = timeit.default_timer()
    
#     with AlignmentFile(bam_file, "rb") as bam, AlignmentFile(out_bam, "wb", header=bam.header) as out_file:
        
#         if chrom not in bam.header.references:
#             print(f"Skipping missing chrom: {chrom}")
#             return
        
#         for read1, read2 in paired_read_gen_stat(bam, read_stats, chrom=chrom):
            
#             # chrom_read_count += 1
#             read_stats.remap_pair += 1
#             og_name = read1.query_name
#             r1_og_seq = read1.query_sequence
#             r1_align_pos = read1.reference_start
#             r2_og_seq = read2.query_sequence
#             r2_align_pos = read2.reference_start
            
#             write_num = 0 # Counter that tracks reads written
            
#             # Get snp df
#             r1_df = r1_het_dict.get(og_name)
#             r2_df = r2_het_dict.get(og_name)
            
            
#             # Og version using a func
#             if r1_df is not None:
#                 r1_het_data = get_read_het_data(r1_df, read1, hap_cols)
                
#                 if r1_het_data is None:
#                     read_stats.discard_unmapped += 1
#                     # SNP overlaps unmapped pos
#                     continue
#                 r1_hap_list = [*make_phased_seqs(r1_het_data[0], *r1_het_data[1])]
            
#             else:
#                 r1_hap_list = [r1_og_seq, r1_og_seq]

            
#             if r2_df is not None:
#                 r2_het_data = get_read_het_data(r2_df, read2, hap_cols)
                
#                 if r2_het_data is None:
#                     read_stats.discard_unmapped += 1
#                     # SNP overlaps unmapped pos
#                     continue
                
#                 r2_hap_list = [*make_phased_seqs(r2_het_data[0], *r2_het_data[1])]

#             else:
#                 r2_hap_list = [r2_og_seq, r2_og_seq]
            
#             # Create pairs to write
#             write_pair_list = [(r1_hap_seq, r2_hap_seq)
#                                for r1_hap_seq, r2_hap_seq in zip(r1_hap_list, r2_hap_list)
#                                if (r1_hap_seq != r1_og_seq) or (r2_hap_seq != r2_og_seq)]
            
#             write_total = len(write_pair_list)
            
#             # Get read pairs
#             for r1_hap_seq, r2_hap_seq in write_pair_list:
#                 write_num += 1
#                 new_read_name = f"{og_name}_WASP_{r1_align_pos}_{r2_align_pos}_{write_num}_{write_total}"
#                 write_read(out_file, read1, r1_hap_seq, new_read_name)
#                 write_read(out_file, read2, r2_hap_seq, new_read_name)
#                 read_stats.write_pair += 1
#                 # chrom_write_count += 1

#         # WOWOW
#         # print(f"{chrom}: Processed {read_stats.remap_pair} pairs and wrote {read_stats.write_pair} new pairs in {timeit.default_timer() - start_chrom:.2f} seconds")
    
#     # Collate and write out fastq now
#     collate_bam = str(Path(out_dir) / f"collate_{chrom}.bam")
#     r1_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r1.fq")
#     r2_out = str(Path(out_dir) / f"swapped_alleles_{chrom}_r2.fq")
    
#     # Do I need to make another file???
#     pysam.collate(out_bam, "-o", collate_bam, catch_stdout=False)
#     pysam.fastq(collate_bam, "-1", r1_out, "-2", r2_out, catch_stdout=False)
#     # print(f"Created fastqs to be remapped in {Path(out_dir) / 'swapped_alleles_{chrom}_r*.fq'}")