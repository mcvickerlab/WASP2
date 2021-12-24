"""
Author: Aaron Ho
Python Version: 3.8
"""


# Default Python package Imports
import time
from collections import Counter

# External package imports
import numpy as np
import pandas as pd
from pandas.arrays import SparseArray
from pysam import VariantFile
from pysam.libcalignmentfile import AlignmentFile


def parse_barcode(bc_series, read):
    """
    Retrieve barcode from read and return grouping

    :param Series bc_series: Barcode group map
    :param PileupRead read: pysam read object
    :return str: Cell type / Cluster
    """
    try:
        barcode = read.alignment.get_tag("CB")
        return bc_series.get(barcode)

    except KeyError:
        return None


def pileup_pos(bam, bc_series, chrom, snp_pos):
    """
    Create pileup column of reads at snp position

    :param AlignmentFile bam: pysam AlignmentFile for bam
    :param str chrom: Chromosome name
    :param int snp_pos: Position of snp in bp
    :return: List of read names and alleles at snp pos
    :rtype: Tuple of (list of str, list of str)
    """
    pile = bam.pileup(chrom, snp_pos-1, snp_pos, truncate=True)

    try:
        pile_col = next(pile)
        return (pile_col.get_query_names(), pile_col.get_query_sequences(),
                [parse_barcode(bc_series, read) for read in pile_col.pileups])

    except StopIteration:
        return None


def count_snp_alleles(bam_file, bc_series, chrom, snp_list, ref_indices, alt_indices):
    """
    Get ref and alt counts of snp's in list

    :param str bam_file: Path to BAM file
    :param str chrom: Chromosome name
    :param snp_list: List of snp tuples
    :type snp_list: list of (int, str, str)
    :return: List of ref count, alt count, other count
    :rtype: List of (int, int, int)
    """
    counted_reads = set()
    allele_counts = []
    
    num_cols = (len(ref_indices) * 2) + 1

    bam = AlignmentFile(bam_file, "rb")

    for snp in snp_list:
        pile_tup = pileup_pos(bam, bc_series, chrom, snp[0])

        if pile_tup is not None:
            read_names, read_alleles, read_groups = pile_tup
            
            count_list = []
            for read_id, allele, group in zip(read_names, read_alleles, read_groups):

                if read_id not in counted_reads:
                    counted_reads.add(read_id)
                    allele = allele.upper()
                    
                    if allele == snp[1]:
                        count_list.append(ref_indices.get(group))
                    elif allele == snp[2]:
                        count_list.append(alt_indices.get(group))
                    else:
                        count_list.append(0)

            if not count_list:
                # allele_counts.append(SparseArray(np.zeros(num_cols), fill_value=0))
                allele_counts.append(np.zeros(num_cols, dtype=np.int32))

            else:
                a_counter = Counter(count_list)
                
                count_array = np.zeros(num_cols)
                count_array[np.fromiter(a_counter.keys(), dtype=np.int32)] = np.fromiter(a_counter.values(), dtype=np.int32)
                
                # allele_counts.append(SparseArray(count_array, fill_value=0))
                allele_counts.append(count_array)

        else:
            # allele_counts.append(SparseArray(np.zeros(num_cols), fill_value=0))
            allele_counts.append(np.zeros(num_cols, dtype=np.int32))

    bam.close()

    return allele_counts


def make_col_data(cell_groups):
    """
    Make column data dynamically from barcode mappings

    :param Series cell_groups: Series containing barcodes as indices, and groupings as items
    :return : list containing list of column names, dict of ref column indices, and dict of alt column indices
    :rtype: Tuple of (list, dict, dict)
    """
    ref_indices = {None: 1}
    alt_indices = {None: 2}
    cols = ["other_count", "noPred_ref", "noPred_alt"]
    
    cell_cols = []
    cell_indices = [i for i in range(3, (len(cell_groups) * 2) + 2, 2)]
    
    for index, cell in zip(cell_indices, cell_groups):
        cell_cols.append(f"{cell}_ref")
        ref_indices[cell] = index
        
        cell_cols.append(f"{cell}_alt")
        alt_indices[cell] = index + 1
    
    cols.extend(cell_cols)
    
    return cols, ref_indices, alt_indices


def make_count_df_sc(bam_file, df, bc_series):
    """
    Make DF containing all intersections and allele counts

    :param str bam_file: Path to BAM file
    :param DataFrame df: Dataframe of intersections, output from
        parse_(intersect/gene)_df()
    :return DataFrame: DataFrame of counts
    """
    count_list = []
    chrom_list = df["chrom"].unique()
    cell_groups = bc_series.unique()
    
    cols, ref_indices, alt_indices = make_col_data(cell_groups)
    skip_chrom = []
    
    total_start = time.time()

    for chrom in chrom_list:
        print(f"Counting Alleles for {chrom}")

        snp_list = df.loc[df["chrom"] == chrom][
            ["pos", "ref", "alt"]].to_records(index=False)

        start = time.time()

        try:
            count_list.extend(count_snp_alleles(bam_file, bc_series, chrom, snp_list, ref_indices, alt_indices))
        except ValueError:
            skip_chrom.append(chrom)
            print(f"Skipping {chrom}: Contig not found\n")
        else:
            print(f"Counted {len(snp_list)} SNP's in {time.time() - start} seconds!\n")

    total_end = time.time()
    print(f"Counted all SNP's in {total_end - total_start} seconds!")

    if skip_chrom:
        df = df.loc[df["chrom"].isin(skip_chrom) == False]

    df[cols] = np.array(count_list, dtype=np.int32)
    df = df.astype({group: "Sparse[int]" for group in cols})
    return df
