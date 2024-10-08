import timeit
from pathlib import Path
from bisect import bisect_left

import polars as pl

from pysam.libcalignmentfile import AlignmentFile

# Helper that does binary search
def find_read_aln_pos(read, pos):
    
    aln_list = read.get_aligned_pairs(True)

    i = bisect_left(aln_list, pos, key=lambda x: x[1])
    
    if i != len(aln_list) and aln_list[i][1] == pos:
        return aln_list[i][0]
    else:
        return None


def make_count_df(bam_file, df):
    """
    Make DF containing all intersections and allele counts

    :param str bam_file: Path to BAM file
    :param DataFrame df: Dataframe of intersections, output from
        parse_(intersect/gene)_df()
    :return DataFrame: DataFrame of counts
    """
    count_list = []

    chrom_list = df.get_column("chrom").unique(
        maintain_order=True)

    total_start = timeit.default_timer()
    
    with AlignmentFile(bam_file, "rb") as bam:
        
        for chrom in chrom_list:
            chrom_df = df.filter(pl.col("chrom") == chrom)
            
            snp_list = chrom_df.select(
                ["pos", "ref", "alt"]).unique(
                subset=["pos"], maintain_order=True).iter_rows()
            
            start = timeit.default_timer()

            try:
                count_list.extend(count_snp_alleles(bam, chrom, snp_list))
            except ValueError:
                print(f"Skipping {chrom}: Contig not found\n")
            else:
                print(f"{chrom}: Counted {chrom_df.height} SNP's in {timeit.default_timer() - start:.2f} seconds!")
                

        total_end = timeit.default_timer()
        print(f"Counted all SNP's in {total_end - total_start:.2f} seconds!")
        
        # Previously used str as chrom instead of cat
        chrom_enum = pl.Enum(df.get_column("chrom").cat.get_categories())
        
        count_df = pl.DataFrame(
            count_list,
            schema={"chrom": chrom_enum,
                    "pos": pl.UInt32,
                    "ref_count": pl.UInt16,
                    "alt_count": pl.UInt16,
                    "other_count": pl.UInt16
                   }
        )
        
        # possibly find better solution
        df = df.with_columns([pl.col("chrom").cast(chrom_enum)]
                             ).join(count_df, on=["chrom", "pos"], how="left")
        
        # df = df.join(count_df, on=["chrom", "pos"], how="left")
    
    return df


def count_snp_alleles(bam, chrom, snp_list):
    """
    Helper function called by...
    make_count_df()
    """
    
    read_set = set()
    allele_counts = []

    for pos, ref, alt in snp_list:

        # read_set = set()
        ref_count, alt_count, other_count = 0, 0, 0

        # Got make sure read is not double counted
        for read in bam.fetch(chrom, pos-1, pos):
            
            # If already counted allele
            if read.query_name in read_set:
                continue
            
            read_set.add(read.query_name)
            
            seq = read.query_sequence
            
            for qpos, refpos in read.get_aligned_pairs(True):
                
                # TODO Update with binary search
                if refpos == pos-1:
                    
                    if seq[qpos] == ref:
                        ref_count+=1
                    elif seq[qpos] == alt:
                        alt_count+=1
                    else:
                        other_count+=1
                    
                    # Found no longer need to loop
                    break
        
        allele_counts.append((chrom, pos, ref_count, alt_count, other_count))
                
    return allele_counts