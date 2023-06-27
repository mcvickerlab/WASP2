
import polars as pl

import pysam
from pysam.libcalignmentfile import AlignmentFile

# Generator for iterating through bam
def paired_read_gen(bam, chrom=None):

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


def paired_read_gen_stat(bam, read_stats, chrom=None):

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


def align_pos_gen(read, align_dict, pos_list):

    yield 0 # yield initial index

    for start, stop in pos_list:
        align_start = align_dict[start]
        
        # for snps, may need to change for indel
        align_stop = align_start + (stop - start)
        
        yield align_start
        yield align_stop
    
    yield len(read.query_sequence)


def get_read_het_data(read_df, read, col_list, max_seqs=None):

    # TODO MULTISAMP AND MAX SEQS
    align_dict = {ref_i: read_i for read_i, ref_i in read.get_aligned_pairs(matches_only=True)}
    pos_list = read_df.select(["start", "stop"]).rows()
    
    try:
        split_pos = [i for i in align_pos_gen(read, align_dict, pos_list)]
        split_seq = [read.query_sequence[start:stop] for start, stop in zip(split_pos[:-1:], split_pos[1:])]
        return split_seq, read_df.select(pl.col(col_list)).get_columns()
    
    except KeyError:
        # remove reads overlap unmapped/gap
        return None


# def get_read_het_data(read_df, read, hap1_col, hap2_col, max_seqs=None):

#     # TODO MULTISAMP AND MAX SEQS
#     align_dict = {ref_i: read_i for read_i, ref_i in read.get_aligned_pairs(matches_only=True)}
#     pos_list = read_df.select(["start", "stop"]).rows()
    
#     try:
#         split_pos = [i for i in align_pos_gen(read, align_dict, pos_list)]
#         split_seq = [read.query_sequence[start:stop] for start, stop in zip(split_pos[:-1:], split_pos[1:])]
#         return split_seq, read_df.get_column(hap1_col), read_df.get_column(hap2_col)
    
#     except KeyError:
#         # remove reads overlap unmapped/gap
#         return None


def make_phased_seqs(split_seq, hap1_alleles, hap2_alleles):
    
    hap1_split = split_seq.copy()
    hap2_split = split_seq.copy()

    hap1_split[1::2] = hap1_alleles
    hap2_split[1::2] = hap2_alleles
    
    return "".join(hap1_split), "".join(hap2_split)


def make_multi_seqs(split_seq, allele_combos):
    
    seq_list = []
    for phased_alleles in allele_combos:
        
        hap_split = split_seq.copy()
        hap_split[1::2] = phased_alleles
        seq_list.append("".join(hap_split))
    
    return seq_list


def write_read(out_bam, read, new_seq, new_name):
    og_qual = read.query_qualities
    read.query_sequence = new_seq
    read.query_name = new_name
    read.query_qualities = og_qual
    out_bam.write(read)