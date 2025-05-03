"""
Author: Aaron Ho
Python Version: 3.8
"""

# Default Python package Imports
import time
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd

# External package imports
import numpy as np
from pysam.libcalignmentfile import AlignmentFile
from pysam import VariantFile


def parse_barcode(bc_series: pd.Series, read: Any) -> Optional[str]:
    """
    Retrieve barcode from a read and return its grouping.

    This function extracts the "CB" tag (cell barcode) from the alignment of a read and then
    uses the provided barcode series to look up the corresponding group (e.g., cell type or cluster).

    Parameters
    ----------
    bc_series : pandas.Series
        Barcode group mapping, with barcodes as keys and group names as values.
    read : pysam.PileupRead
        A pysam read object from which to extract the barcode.

    Returns
    -------
    str or None
        The group (e.g., cell type or cluster) corresponding to the read's barcode, or None if not found.

    Examples
    --------
    >>> group = parse_barcode(bc_series, read)
    """
    try:
        barcode = read.alignment.get_tag("CB")
        return bc_series.get(barcode)
    except KeyError:
        return None


def pileup_pos(bam: AlignmentFile, bc_series: pd.Series, chrom: str, snp_pos: int) -> Optional[Tuple[List[str], List[str], List[Optional[str]]]]:
    """
    Create a pileup column of reads at a SNP position and retrieve barcode groupings.

    This function generates a pileup from the BAM file for the specified chromosome and SNP position.
    It returns a tuple containing the read names, query sequences, and the corresponding barcode groupings
    (obtained via the `parse_barcode` function) for each read in the pileup.

    Parameters
    ----------
    bam : AlignmentFile
        A pysam AlignmentFile instance for the BAM file.
    bc_series : pandas.Series
        A mapping of barcodes to group names.
    chrom : str
        Chromosome name.
    snp_pos : int
        Position of the SNP in base pairs.

    Returns
    -------
    tuple of (list of str, list of str, list of str) or None
        A tuple containing:
            - A list of read names.
            - A list of query sequences at the SNP position.
            - A list of barcode groupings for each read.

        Returns None if no pileup is generated.

    Examples
    --------
    >>> result = pileup_pos(bam, bc_series, "chr1", 123456)
    """
    pile = bam.pileup(chrom, snp_pos - 1, snp_pos, truncate=True)

    try:
        pile_col = next(pile)
        return (
            pile_col.get_query_names(),
            pile_col.get_query_sequences(),
            [parse_barcode(bc_series, read) for read in pile_col.pileups]
        )
    except StopIteration:
        return None


def count_snp_alleles(bam_file: str, bc_series: pd.Series, chrom: str, snp_list: List[Tuple[int, str, str]], 
                      ref_indices: Dict[Any, int], alt_indices: Dict[Any, int]) -> List[np.ndarray]:
    """
    Get reference and alternate allele counts for a list of SNPs.

    This function opens the specified BAM file and, for each SNP in the provided list, creates a pileup
    at the SNP position to count the occurrence of reference and alternate alleles. The allele is tallied
    based on barcode groupings using the provided reference and alternate indices. If a read's allele does
    not match either the reference or alternate, it is counted as "other".

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    bc_series : pandas.Series
        Barcode group mapping to group names.
    chrom : str
        Chromosome name.
    snp_list : list of tuple
        List of SNP tuples, where each tuple contains:
            (position (int), reference allele (str), alternate allele (str)).
    ref_indices : dict
        Dictionary mapping group names to reference column indices.
    alt_indices : dict
        Dictionary mapping group names to alternate column indices.

    Returns
    -------
    list of numpy.ndarray
        A list where each element is a numpy array of allele counts. The number of columns is
        computed as (len(ref_indices) * 2) + 1. Typically, column 0 represents "other" counts,
        while the remaining columns represent counts for reference and alternate alleles for each group.

    Examples
    --------
    >>> counts = count_snp_alleles("example.bam", bc_series, "chr1", snp_list, ref_indices, alt_indices)
    """
    counted_reads = set()
    allele_counts: List[np.ndarray] = []
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
                allele_counts.append(np.zeros(num_cols, dtype=np.int32))
            else:
                a_counter = Counter(count_list)
                count_array = np.zeros(num_cols, dtype=np.int32)
                count_array[np.fromiter(a_counter.keys(), dtype=np.int32)] = np.fromiter(a_counter.values(), dtype=np.int32)
                allele_counts.append(count_array)
        else:
            allele_counts.append(np.zeros(num_cols, dtype=np.int32))

    bam.close()

    return allele_counts


def make_col_data(cell_groups: pd.Series) -> Tuple[List[str], Dict[Any, int], Dict[Any, int]]:
    """
    Create column metadata dynamically from barcode mappings.

    This function generates a list of column names and corresponding dictionaries for reference
    and alternate allele indices based on the unique cell groups provided.

    Parameters
    ----------
    cell_groups : pandas.Series
        Series containing barcodes as indices and group names as values.

    Returns
    -------
    tuple
        A tuple containing:
            - A list of column names.
            - A dictionary mapping cell groups to their reference column indices.
            - A dictionary mapping cell groups to their alternate column indices.

    Examples
    --------
    >>> cols, ref_indices, alt_indices = make_col_data(cell_groups)
    """
    ref_indices: Dict[Any, int] = {None: 1}
    alt_indices: Dict[Any, int] = {None: 2}
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


def make_count_df_sc(bam_file: str, df: pd.DataFrame, bc_series: pd.Series) -> pd.DataFrame:
    """
    Create a DataFrame with allele counts for single-cell data.

    This function iterates over each unique chromosome in the provided DataFrame and uses the specified
    BAM file along with a barcode series to count allele occurrences for each SNP. It dynamically creates
    column names based on barcode mappings and appends the allele counts to the DataFrame. Any chromosomes
    not found in the BAM file are skipped.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    df : pandas.DataFrame
        DataFrame of intersections (typically output from parse_intersect_df() or parse_gene_df())
        that includes at least the columns "chrom", "pos", "ref", and "alt".
    bc_series : pandas.Series
        Series mapping barcodes to their respective cell groups.

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with additional columns for allele counts. The columns include:
            - Reference counts.
            - Alternate counts.
            - Other counts.

        The resulting DataFrame is cast to use Sparse integer arrays for the allele count columns.

    Examples
    --------
    >>> df_counts = make_count_df_sc("example.bam", df_intersections, bc_series)
    """
    count_list = []
    chrom_list = df["chrom"].unique()
    cell_groups = bc_series.unique()
    
    cols, ref_indices, alt_indices = make_col_data(cell_groups)
    skip_chrom = []
    
    total_start = time.time()

    for chrom in chrom_list:
        print(f"Counting Alleles for {chrom}")

        snp_list = df.loc[df["chrom"] == chrom][["pos", "ref", "alt"]].to_records(index=False)

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
