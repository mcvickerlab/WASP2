from typing import Optional, Generator, Tuple, Dict, List, Any
import numpy as np

import polars as pl

import pysam
from pysam import AlignmentFile, AlignedSegment

# Generator for iterating through bam
def paired_read_gen(
    bam: AlignmentFile,
    chrom: Optional[str] = None
) -> Generator[Tuple[AlignedSegment, AlignedSegment], None, None]:

    read_dict = {}
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


def paired_read_gen_stat(
    bam: AlignmentFile,
    read_stats: Any,
    chrom: Optional[str] = None
) -> Generator[Tuple[AlignedSegment, AlignedSegment], None, None]:

    read_dict = {}
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


def align_pos_gen(
    read: AlignedSegment,
    align_dict: Dict[int, int],
    pos_list: List[Tuple[int, int]]
) -> Generator[int, None, None]:

    yield 0 # yield initial index

    for start, stop in pos_list:
        align_start = align_dict[start]
        
        # for snps, may need to change for indel
        align_stop = align_start + (stop - start)
        
        yield align_start
        yield align_stop
    
    yield len(read.query_sequence)


def _build_ref2read_maps(read: AlignedSegment) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Build reference position to read position mappings for indel support.

    Args:
        read: pysam AlignedSegment

    Returns:
        Tuple of (ref2q_left, ref2q_right) dictionaries mapping reference positions
        to read query positions. For deletions (ref pos with no read pos), uses
        nearest left/right query positions.
    """
    # Get all aligned pairs including gaps (matches_only=False)
    # Returns list of (query_pos, ref_pos) tuples, with None for gaps
    pairs = read.get_aligned_pairs(matches_only=False)

    ref2q_left = {}   # Maps ref pos to nearest left query pos
    ref2q_right = {}  # Maps ref pos to nearest right query pos

    last_query_pos = None

    # Forward pass: build left mapping
    for query_pos, ref_pos in pairs:
        if ref_pos is not None:
            if query_pos is not None:
                ref2q_left[ref_pos] = query_pos
                last_query_pos = query_pos
            else:
                # Deletion: use last known query position
                if last_query_pos is not None:
                    ref2q_left[ref_pos] = last_query_pos

    # Backward pass: build right mapping
    last_query_pos = None
    for query_pos, ref_pos in reversed(pairs):
        if ref_pos is not None:
            if query_pos is not None:
                ref2q_right[ref_pos] = query_pos
                last_query_pos = query_pos
            else:
                # Deletion: use next known query position
                if last_query_pos is not None:
                    ref2q_right[ref_pos] = last_query_pos

    return ref2q_left, ref2q_right


def get_read_het_data(
    read_df: pl.DataFrame,
    read: AlignedSegment,
    col_list: List[str],
    max_seqs: Optional[int] = None,
    include_indels: bool = False,
    insert_qual: int = 30
) -> Optional[Tuple[List[str], List[str], List[pl.Series]]]:
    """Extract heterozygous variant data from read with indel support.

    Args:
        read_df: DataFrame with variant positions and alleles
        read: pysam AlignedSegment
        col_list: List of column names containing alleles
        max_seqs: Maximum number of alternate sequences (unused currently)
        include_indels: Whether to use indel-aware position mapping
        insert_qual: Quality score for inserted bases (Phred scale)

    Returns:
        Tuple of (split_seq, split_qual, allele_series) or None if mapping fails
        split_seq: List of sequence segments between variants
        split_qual: List of quality score segments
        allele_series: List of polars Series with allele data
    """
    pos_list = read_df.select(["start", "stop"]).rows()

    try:
        if include_indels:
            # Use indel-aware mapping
            ref2q_left, ref2q_right = _build_ref2read_maps(read)

            split_pos = [0]  # Start with query position 0
            split_qual_pos = [0]

            for start, stop in pos_list:
                # Use left mapping for variant start, right mapping for variant end
                if start not in ref2q_left or stop not in ref2q_right:
                    # Variant overlaps unmapped region
                    return None

                query_start = ref2q_left[start]
                query_stop = ref2q_right[stop]

                split_pos.append(query_start)
                split_pos.append(query_stop)
                split_qual_pos.append(query_start)
                split_qual_pos.append(query_stop)

            split_pos.append(len(read.query_sequence))
            split_qual_pos.append(len(read.query_qualities))

        else:
            # Original SNP-only logic (backward compatible)
            align_dict = {ref_i: read_i for read_i, ref_i in read.get_aligned_pairs(matches_only=True)}
            split_pos = [i for i in align_pos_gen(read, align_dict, pos_list)]
            split_qual_pos = split_pos.copy()

        # Extract sequence and quality segments
        split_seq = [read.query_sequence[start:stop] for start, stop in zip(split_pos[:-1], split_pos[1:])]
        split_qual = [read.query_qualities[start:stop] for start, stop in zip(split_qual_pos[:-1], split_qual_pos[1:])]

        return split_seq, split_qual, read_df.select(pl.col(col_list)).get_columns()

    except KeyError:
        # remove reads overlap unmapped/gap
        return None



def _fill_insertion_quals(insert_len: int, left_qual: np.ndarray, right_qual: np.ndarray, insert_qual: int = 30) -> np.ndarray:
    """Generate quality scores for inserted bases.

    Args:
        insert_len: Number of inserted bases needing quality scores
        left_qual: Quality scores from left flanking region
        right_qual: Quality scores from right flanking region
        insert_qual: Default quality score if flanks unavailable

    Returns:
        Numpy array of quality scores for inserted bases
    """
    if len(left_qual) == 0 and len(right_qual) == 0:
        # No flanking quality data, use constant
        return np.full(insert_len, insert_qual, dtype=np.uint8)

    # Average flanking qualities
    flank_quals = np.concatenate([left_qual, right_qual])
    mean_qual = int(np.mean(flank_quals))
    return np.full(insert_len, mean_qual, dtype=np.uint8)


def make_phased_seqs(split_seq: List[str], hap1_alleles: Any, hap2_alleles: Any) -> Tuple[str, str]:
    """Create phased sequences by swapping alleles (SNP-only version).

    Args:
        split_seq: List of sequence segments
        hap1_alleles: Haplotype 1 alleles
        hap2_alleles: Haplotype 2 alleles

    Returns:
        Tuple of (hap1_seq, hap2_seq) strings
    """
    hap1_split = split_seq.copy()
    hap2_split = split_seq.copy()

    hap1_split[1::2] = hap1_alleles
    hap2_split[1::2] = hap2_alleles

    return "".join(hap1_split), "".join(hap2_split)


def make_phased_seqs_with_qual(
    split_seq: List[str],
    split_qual: List[np.ndarray],
    hap1_alleles: Any,
    hap2_alleles: Any,
    insert_qual: int = 30
) -> Tuple[Tuple[str, np.ndarray], Tuple[str, np.ndarray]]:
    """Create phased sequences with quality scores (indel-aware version).

    Args:
        split_seq: List of sequence segments
        split_qual: List of quality score arrays
        hap1_alleles: Haplotype 1 alleles
        hap2_alleles: Haplotype 2 alleles
        insert_qual: Quality score for inserted bases

    Returns:
        Tuple of ((hap1_seq, hap1_qual), (hap2_seq, hap2_qual))
    """
    hap1_seq_parts = []
    hap1_qual_parts = []
    hap2_seq_parts = []
    hap2_qual_parts = []

    for i, (seq_part, qual_part) in enumerate(zip(split_seq, split_qual)):
        if i % 2 == 0:
            # Non-variant segment - same for both haplotypes
            hap1_seq_parts.append(seq_part)
            hap1_qual_parts.append(qual_part)
            hap2_seq_parts.append(seq_part)
            hap2_qual_parts.append(qual_part)
        else:
            # Variant segment - swap alleles
            idx = i // 2
            hap1_allele = hap1_alleles[idx]
            hap2_allele = hap2_alleles[idx]

            hap1_seq_parts.append(hap1_allele)
            hap2_seq_parts.append(hap2_allele)

            # Handle quality scores for insertions/deletions
            orig_len = len(seq_part)
            hap1_len = len(hap1_allele)
            hap2_len = len(hap2_allele)

            # Get flanking quality scores for insertion quality inference
            left_qual = split_qual[i-1] if i > 0 else np.array([], dtype=np.uint8)
            right_qual = split_qual[i+1] if i < len(split_qual) - 1 else np.array([], dtype=np.uint8)

            # Haplotype 1 quality handling
            if hap1_len == orig_len:
                # Same length - use original qualities
                hap1_qual_parts.append(qual_part)
            elif hap1_len < orig_len:
                # Deletion - truncate qualities
                hap1_qual_parts.append(qual_part[:hap1_len])
            else:
                # Insertion - fill extra qualities
                extra_len = hap1_len - orig_len
                extra_quals = _fill_insertion_quals(extra_len, left_qual, right_qual, insert_qual)
                hap1_qual_parts.append(np.concatenate([qual_part, extra_quals]))

            # Haplotype 2 quality handling
            if hap2_len == orig_len:
                hap2_qual_parts.append(qual_part)
            elif hap2_len < orig_len:
                hap2_qual_parts.append(qual_part[:hap2_len])
            else:
                extra_len = hap2_len - orig_len
                extra_quals = _fill_insertion_quals(extra_len, left_qual, right_qual, insert_qual)
                hap2_qual_parts.append(np.concatenate([qual_part, extra_quals]))

    hap1_seq = "".join(hap1_seq_parts)
    hap2_seq = "".join(hap2_seq_parts)
    hap1_qual = np.concatenate(hap1_qual_parts)
    hap2_qual = np.concatenate(hap2_qual_parts)

    return (hap1_seq, hap1_qual), (hap2_seq, hap2_qual)


def make_multi_seqs(split_seq: List[str], allele_combos: Any) -> List[str]:
    """Create multiple sequences for multi-sample analysis (SNP-only version).

    Args:
        split_seq: List of sequence segments
        allele_combos: List of allele combinations across samples

    Returns:
        List of sequence strings, one per unique haplotype
    """
    seq_list = []
    for phased_alleles in allele_combos:

        hap_split = split_seq.copy()
        hap_split[1::2] = phased_alleles
        seq_list.append("".join(hap_split))

    return seq_list


def make_multi_seqs_with_qual(
    split_seq: List[str],
    split_qual: List[np.ndarray],
    allele_combos: Any,
    insert_qual: int = 30
) -> List[Tuple[str, np.ndarray]]:
    """Create multiple sequences with quality scores for multi-sample indel support.

    Args:
        split_seq: List of sequence segments
        split_qual: List of quality score arrays
        allele_combos: List of allele combinations across samples
        insert_qual: Quality score for inserted bases

    Returns:
        List of (sequence, quality) tuples, one per unique haplotype
    """
    result_list = []

    for phased_alleles in allele_combos:
        seq_parts = []
        qual_parts = []

        for i, (seq_part, qual_part) in enumerate(zip(split_seq, split_qual)):
            if i % 2 == 0:
                # Non-variant segment - use as is
                seq_parts.append(seq_part)
                qual_parts.append(qual_part)
            else:
                # Variant segment - use allele from this haplotype
                idx = i // 2
                allele = phased_alleles[idx]
                seq_parts.append(allele)

                # Handle quality scores for length differences
                orig_len = len(seq_part)
                allele_len = len(allele)

                # Get flanking qualities
                left_qual = split_qual[i-1] if i > 0 else np.array([], dtype=np.uint8)
                right_qual = split_qual[i+1] if i < len(split_qual) - 1 else np.array([], dtype=np.uint8)

                if allele_len == orig_len:
                    # Same length - use original qualities
                    qual_parts.append(qual_part)
                elif allele_len < orig_len:
                    # Deletion - truncate qualities
                    qual_parts.append(qual_part[:allele_len])
                else:
                    # Insertion - fill extra qualities
                    extra_len = allele_len - orig_len
                    extra_quals = _fill_insertion_quals(extra_len, left_qual, right_qual, insert_qual)
                    qual_parts.append(np.concatenate([qual_part, extra_quals]))

        hap_seq = "".join(seq_parts)
        hap_qual = np.concatenate(qual_parts)
        result_list.append((hap_seq, hap_qual))

    return result_list


def write_read(out_bam: AlignmentFile, read: AlignedSegment, new_seq: str, new_name: str, new_qual: Optional[np.ndarray] = None) -> None:
    """Write a modified read to output BAM.

    Args:
        out_bam: Output BAM file
        read: Original read
        new_seq: New sequence
        new_name: New read name
        new_qual: Optional new quality scores (for indels)
    """
    if new_qual is None:
        # SNP mode - preserve original qualities (sequence length unchanged)
        og_qual = read.query_qualities
        read.query_sequence = new_seq
        read.query_name = new_name
        read.query_qualities = og_qual
    else:
        # Indel mode - use provided qualities
        # CIGAR must match sequence length, update if length changed
        old_len = read.query_length
        new_len = len(new_seq)
        if old_len != new_len:
            # Sequence length changed due to indel, update CIGAR to simple match
            # These reads will be realigned anyway during remapping
            read.cigartuples = [(0, new_len)]  # 0 = MATCH operation
        read.query_sequence = new_seq
        read.query_name = new_name
        read.query_qualities = new_qual
    out_bam.write(read)