"""
Author: Aaron Ho
Python Version: 3.8
"""

# Default Python package Imports
import logging
import time
from collections import Counter

# External package imports
from pysam.libcalignmentfile import AlignmentFile

logger = logging.getLogger(__name__)


def pileup_pos(bam, chrom, snp_pos):
    """
    Create pileup column of reads at snp position

    :param AlignmentFile bam: pysam AlignmentFile for bam
    :param str chrom: Chromosome name
    :param int snp_pos: Position of snp in bp
    :return: List of read names and alleles at snp pos
    :rtype: Tuple of (list of str, list of str)
    """
    pile = bam.pileup(chrom, snp_pos - 1, snp_pos, truncate=True)

    try:
        pile_col = next(pile)
        return pile_col.get_query_names(), pile_col.get_query_sequences()

    except StopIteration:
        logger.debug("No pileup data at %s:%d", chrom, snp_pos)
        return None


def count_snp_alleles(bam_file, chrom, snp_list):
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


def make_count_df(bam_file, df):
    """
    Make DF containing all intersections and allele counts

    :param str bam_file: Path to BAM file
    :param DataFrame df: Dataframe of intersections, output from
        parse_(intersect/gene)_df()
    :return DataFrame: DataFrame of counts
    """
    count_list = []
    chrom_list = df["chrom"].unique()
    skip_chrom = []

    total_start = time.time()

    for chrom in chrom_list:
        logger.info("Counting alleles for %s", chrom)

        snp_list = df.loc[df["chrom"] == chrom][["pos", "ref", "alt"]].to_records(index=False)

        start = time.time()

        try:
            count_list.extend(count_snp_alleles(bam_file, chrom, snp_list))
        except ValueError:
            skip_chrom.append(chrom)
            logger.warning("Skipping %s: contig not found", chrom)
        else:
            logger.info("Counted %d SNPs in %.2f seconds", len(snp_list), time.time() - start)

    total_end = time.time()
    logger.info("Counted all SNPs in %.2f seconds", total_end - total_start)

    if skip_chrom:
        df = df.loc[not df["chrom"].isin(skip_chrom)]

    df[["ref_count", "alt_count", "other_count"]] = count_list
    return df
