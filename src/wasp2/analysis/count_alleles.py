"""
Author: Aaron Ho
Python Version: 3.8
"""

# Default Python package Imports
import time
from collections import Counter
from typing import List, Optional, Tuple
import pandas as pd

# External package imports
from pysam.libcalignmentfile import AlignmentFile


def pileup_pos(bam: AlignmentFile, chrom: str, snp_pos: int) -> Optional[Tuple[List[str], List[str]]]:
    """
    Create a pileup column of reads at a SNP position.

    This function generates a pileup from the given BAM file for a specified chromosome
    and SNP position, then returns a tuple containing the read names and corresponding
    query sequences at that position.

    Parameters
    ----------
    bam : AlignmentFile
        A pysam AlignmentFile instance representing the BAM file.
    chrom : str
        The name of the chromosome.
    snp_pos : int
        The base pair position of the SNP.

    Returns
    -------
    tuple of (list of str, list of str) or None
        A tuple containing:
            - A list of read names.
            - A list of query sequences at the SNP position.

        Returns None if no pileup is available.

    Examples
    --------
    >>> bam = AlignmentFile("example.bam", "rb")
    >>> names, sequences = pileup_pos(bam, "chr1", 123456)
    """
    pile = bam.pileup(chrom, snp_pos - 1, snp_pos, truncate=True)

    try:
        pile_col = next(pile)
        return pile_col.get_query_names(), pile_col.get_query_sequences()
    except StopIteration:
        return None


def count_snp_alleles(bam_file: str, chrom: str, snp_list: List[Tuple[int, str, str]]) -> List[Tuple[int, int, int]]:
    """
    Get reference and alternate allele counts for SNPs in a list.

    This function opens a BAM file and, for each SNP in the provided list, creates a pileup
    at the SNP position to count the occurrence of reference and alternate alleles.
    It returns a list of allele counts as tuples (ref_count, alt_count, other_count).

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    chrom : str
        The name of the chromosome.
    snp_list : list of tuple
        A list of SNP tuples, where each tuple contains:
            (position (int), reference allele (str), alternate allele (str)).

    Returns
    -------
    list of tuple
        A list where each element is a tuple of counts:
            (ref_count, alt_count, other_count).

    Examples
    --------
    >>> snps = [(123456, "A", "T"), (234567, "C", "G")]
    >>> counts = count_snp_alleles("example.bam", "chr1", snps)
    """
    counted_reads = set()
    allele_counts: List[Tuple[int, int, int]] = []

    bam = AlignmentFile(bam_file, "rb")

    for snp in snp_list:
        pile_tup = pileup_pos(bam, chrom, snp[0])

        if pile_tup is not None:
            read_names, read_alleles = pile_tup
            count_list = []

            for read_id, allele in zip(read_names, read_alleles):
                if read_id not in counted_reads:
                    counted_reads.add(read_id)
                    count_list.append(allele.upper())

            if not count_list:
                allele_counts.append((0, 0, 0))
            else:
                a_counter = Counter(count_list)
                total_count = sum(a_counter.values())

                ref_count = a_counter.get(snp[1], 0)
                alt_count = a_counter.get(snp[2], 0)

                allele_counts.append((ref_count, alt_count, total_count - ref_count - alt_count))
        else:
            allele_counts.append((0, 0, 0))

    bam.close()

    return allele_counts


def make_count_df(bam_file: str, df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a DataFrame containing allele counts for SNP intersections.

    This function iterates over each unique chromosome in the provided DataFrame, counts the alleles
    at SNP positions using the specified BAM file, and appends the counts as new columns to the DataFrame.
    If certain chromosomes are not found in the BAM file, they are skipped.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    df : pandas.DataFrame
        A DataFrame of intersections (typically output from parse_intersect_df() or parse_gene_df())
        with at least the columns "chrom", "pos", "ref", and "alt".

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with additional columns:
            - "ref_count": Count of reference alleles.
            - "alt_count": Count of alternate alleles.
            - "other_count": Count of alleles that are neither ref nor alt.

    Examples
    --------
    >>> df_counts = make_count_df("example.bam", df_intersections)
    """
    count_list = []
    chrom_list = df["chrom"].unique()
    skip_chrom = []

    total_start = time.time()

    for chrom in chrom_list:
        print(f"Counting Alleles for {chrom}")

        snp_list = df.loc[df["chrom"] == chrom][
            ["pos", "ref", "alt"]].to_records(index=False)

        start = time.time()

        try:
            count_list.extend(count_snp_alleles(bam_file, chrom, snp_list))
        except ValueError:
            skip_chrom.append(chrom)
            print(f"Skipping {chrom}: Contig not found\n")
        else:
            print(f"Counted {len(snp_list)} SNP's in {time.time() - start} seconds!\n")

    total_end = time.time()
    print(f"Counted all SNP's in {total_end - total_start} seconds!")

    if skip_chrom:
        df = df.loc[df["chrom"].isin(skip_chrom) == False]

    df[["ref_count", "alt_count", "other_count"]] = count_list
    return df
