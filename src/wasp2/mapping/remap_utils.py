from typing import List, Optional, Tuple, Iterable, Dict, Any

import polars as pl

import pysam
from pysam.libcalignmentfile import AlignmentFile

# Generator for iterating through BAM
def paired_read_gen(bam: AlignmentFile, chrom: Optional[str] = None) -> Iterable[Tuple[pysam.AlignedSegment, pysam.AlignedSegment]]:
    """
    Generator that yields paired reads from a BAM file.

    Iterates through reads fetched from the BAM file (optionally for a specific chromosome)
    and yields a tuple of read pairs. Only proper pairs (i.e. not secondary or supplementary)
    are returned.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        An open BAM file.
    chrom : str, optional
        Chromosome to filter reads by. If None, all chromosomes are processed.

    Yields
    ------
    tuple
        A tuple of two pysam.AlignedSegment objects representing a paired read.
    """
    read_dict: Dict[str, pysam.AlignedSegment] = {}
    for read in bam.fetch(chrom):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue

        if read.query_name not in read_dict:
            read_dict[read.query_name] = read
            continue

        if read.is_read1:
            yield read, read_dict.pop(read.query_name)
        else:
            yield read_dict.pop(read.query_name), read


def paired_read_gen_stat(bam: AlignmentFile, read_stats: Any, chrom: Optional[str] = None) -> Iterable[Tuple[pysam.AlignedSegment, pysam.AlignedSegment]]:
    """
    Generator that yields paired reads from a BAM file and updates read statistics.

    Iterates through reads from the BAM file (optionally for a specific chromosome) and yields
    paired reads. It also updates the provided read_stats object with counts of discarded reads
    based on improper pairing, secondary, or supplementary flags.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        An open BAM file.
    read_stats : Any
        An object to track statistics about discarded and processed reads.
    chrom : str, optional
        Chromosome to filter reads by. If None, all chromosomes are processed.

    Yields
    ------
    tuple
        A tuple of two pysam.AlignedSegment objects representing a paired read.

    Notes
    -----
    After processing, the number of missing pairs is added to read_stats.discard_missing_pair.
    """
    read_dict: Dict[str, pysam.AlignedSegment] = {}
    discard_set = set()
    
    # DO I need multiple iterators???
    for read in bam.fetch(chrom, multiple_iterators=False):
        if not read.is_proper_pair:
            discard_set.add(read.query_name)
            read_stats.discard_improper_pair += 1
            continue
        elif read.is_secondary:
            discard_set.add(read.query_name)
            read_stats.discard_secondary += 1
            continue
        elif read.is_supplementary:
            discard_set.add(read.query_name)
            read_stats.discard_supplementary += 1
            continue

        if read.query_name not in read_dict:
            read_dict[read.query_name] = read
            continue

        if read.is_read1:
            yield read, read_dict.pop(read.query_name)
        else:
            yield read_dict.pop(read.query_name), read
    
    # Process missing pairs
    read_stats.discard_missing_pair += len(set(read_dict.keys()) - discard_set)


def align_pos_gen(read: pysam.AlignedSegment, align_dict: Dict[int, int], pos_list: Iterable[Tuple[int, int]]) -> Iterable[int]:
    """
    Generate alignment positions for given intervals within a read.

    Yields the initial index (0), then for each (start, stop) interval provided in pos_list,
    yields the aligned start and calculated aligned stop positions based on the alignment mapping,
    and finally yields the length of the read's query sequence.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to process.
    align_dict : dict
        Dictionary mapping reference positions to query positions from the read's aligned pairs.
    pos_list : iterable
        Iterable of (start, stop) tuples representing intervals in the reference.

    Yields
    ------
    int
        An integer position within the read's query sequence.
    """
    yield 0  # yield initial index

    for start, stop in pos_list:
        align_start = align_dict[start]
        # for SNPs, may need to change for indels
        align_stop = align_start + (stop - start)
        yield align_start
        yield align_stop

    yield len(read.query_sequence)


def get_read_het_data(read_df: pl.DataFrame, read: pysam.AlignedSegment, col_list: List[str], max_seqs: Optional[int] = None) -> Optional[Tuple[List[str], List[pl.Series]]]:
    """
    Extract heterozygous read data from a DataFrame for a given read.

    This function uses the aligned pairs from the read to construct a list of split sequences
    based on positions specified in read_df, and retrieves corresponding genotype data from the DataFrame
    columns specified in col_list.

    Parameters
    ----------
    read_df : polars.DataFrame
        A DataFrame containing columns "start" and "stop", and additional genotype information.
    read : pysam.AlignedSegment
        The read from which to extract the sequence.
    col_list : list
        List of column names in read_df to retrieve genotype data from.
    max_seqs : int, optional
        Maximum number of sequences to process (not currently used).

    Returns
    -------
    tuple or None
        A tuple containing:
            - A list of split sequences derived from the read.
            - A list of columns (as Polars Series) corresponding to col_list.
            
        Returns None if a KeyError occurs (e.g., if the read overlaps an unmapped or gapped region).
    """
    # TODO MULTISAMP AND MAX SEQS
    align_dict: Dict[int, int] = {ref_i: read_i for read_i, ref_i in read.get_aligned_pairs(matches_only=True)}
    pos_list = read_df.select(["start", "stop"]).rows()
    
    try:
        split_pos = [i for i in align_pos_gen(read, align_dict, pos_list)]
        split_seq = [read.query_sequence[start:stop] for start, stop in zip(split_pos[:-1], split_pos[1:])]
        return split_seq, read_df.select(pl.col(col_list)).get_columns()
    
    except KeyError:
        # remove reads overlap unmapped/gap
        return None


# def get_read_het_data(read_df, read, hap1_col, hap2_col, max_seqs=None):
#     # TODO MULTISAMP AND MAX SEQS
#     align_dict = {ref_i: read_i for read_i, ref_i in read.get_aligned_pairs(matches_only=True)}
#     pos_list = read_df.select(["start", "stop"]).rows()
#     
#     try:
#         split_pos = [i for i in align_pos_gen(read, align_dict, pos_list)]
#         split_seq = [read.query_sequence[start:stop] for start, stop in zip(split_pos[:-1:], split_pos[1:])]
#         return split_seq, read_df.get_column(hap1_col), read_df.get_column(hap2_col)
#     
#     except KeyError:
#         # remove reads overlap unmapped/gap
#         return None


def make_phased_seqs(split_seq: List[str], hap1_alleles: List[str], hap2_alleles: List[str]) -> Tuple[str, str]:
    """
    Construct two phased sequences by replacing alternating positions in the split sequence.

    This function creates two sequences by replacing every other segment in the split sequence with
    the provided haplotype alleles (hap1_alleles and hap2_alleles).

    Parameters
    ----------
    split_seq : list of str
        List of sequence segments from the original read.
    hap1_alleles : list of str
        Alleles to insert into the first phased sequence.
    hap2_alleles : list of str
        Alleles to insert into the second phased sequence.

    Returns
    -------
    tuple of str
        A tuple containing the two phased sequences as strings.
    """
    hap1_split = split_seq.copy()
    hap2_split = split_seq.copy()

    hap1_split[1::2] = hap1_alleles
    hap2_split[1::2] = hap2_alleles
    
    return "".join(hap1_split), "".join(hap2_split)


def make_multi_seqs(split_seq: List[str], allele_combos: Iterable[List[str]]) -> List[str]:
    """
    Construct multiple sequences based on different allele combinations.

    For each set of alleles in allele_combos, replaces the alternating segments in the split sequence
    with the given alleles and returns a list of the resulting sequences.

    Parameters
    ----------
    split_seq : list of str
        List of sequence segments from the original read.
    allele_combos : iterable
        An iterable of allele combinations. Each element is used to replace alternating segments.

    Returns
    -------
    list of str
        A list of sequences generated from the given allele combinations.
    """
    seq_list: List[str] = []
    for phased_alleles in allele_combos:
        hap_split = split_seq.copy()
        hap_split[1::2] = phased_alleles
        seq_list.append("".join(hap_split))
    return seq_list


def write_read(out_bam: AlignmentFile, read: pysam.AlignedSegment, new_seq: str, new_name: str) -> None:
    """
    Write a modified read to an output BAM file.

    This function replaces the read's sequence with new_seq and its query name with new_name,
    while preserving the original quality scores, and then writes the modified read to the
    output BAM file.

    Parameters
    ----------
    out_bam : pysam.AlignmentFile
        An open BAM file in write mode.
    read : pysam.AlignedSegment
        The original read.
    new_seq : str
        The new sequence to assign to the read.
    new_name : str
        The new query name for the read.

    Returns
    -------
    None
    """
    og_qual = read.query_qualities
    read.query_sequence = new_seq
    read.query_name = new_name
    read.query_qualities = og_qual
    out_bam.write(read)
