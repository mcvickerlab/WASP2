import timeit
from pathlib import Path
from bisect import bisect_left
from typing import Any, Callable, Iterable, List, Optional, Tuple

import polars as pl

from pysam.libcalignmentfile import AlignmentFile


def find_read_aln_pos(read: Any, pos: int) -> Optional[int]:
    """
    Find the query position in a read that aligns to a given genomic position using binary search.

    This function retrieves aligned pairs from the read and uses a binary search (via :func:`bisect_left`)
    to find the query position corresponding to the given genomic position.

    Parameters
    ----------
    read : pysam.AlignedSegment
        A read object from which to obtain aligned pairs.
    pos : int
        The genomic position to search for.

    Returns
    -------
    int or None
        The query position corresponding to the genomic position if found; otherwise, None.
    """
    aln_list = read.get_aligned_pairs(True)

    i = bisect_left(aln_list, pos, key=lambda x: x[1])
    
    if i != len(aln_list) and aln_list[i][1] == pos:
        return aln_list[i][0]
    else:
        return None


def make_count_df(bam_file: str, df: pl.DataFrame) -> pl.DataFrame:
    """
    Create a DataFrame containing intersections and allele counts.

    This function processes a BAM file and a Polars DataFrame of intersections (e.g., output from
    a parsing function such as :func:`parse_(intersect/gene)_df`). It computes allele counts for each SNP,
    assembles the counts into a new DataFrame, and then joins it with the input DataFrame.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    df : polars.DataFrame
        DataFrame of intersections containing at least a "chrom" column and SNP information.

    Returns
    -------
    polars.DataFrame
        A DataFrame containing the original intersection data joined with computed allele counts.
    """
    count_list: List[Tuple[str, int, int, int, int]] = []

    chrom_list: List[str] = df.get_column("chrom").unique(maintain_order=True)

    total_start = timeit.default_timer()
    
    with AlignmentFile(bam_file, "rb") as bam:
        for chrom in chrom_list:
            chrom_df = df.filter(pl.col("chrom") == chrom)
            
            snp_list = chrom_df.select(["pos", "ref", "alt"]).unique(
                subset=["pos"], maintain_order=True
            ).iter_rows()
            
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
            schema={
                "chrom": chrom_enum,
                "pos": pl.UInt32,
                "ref_count": pl.UInt16,
                "alt_count": pl.UInt16,
                "other_count": pl.UInt16,
            }
        )
        
        # possibly find better solution
        df = df.with_columns([pl.col("chrom").cast(chrom_enum)]).join(count_df, on=["chrom", "pos"], how="left")
        
        # df = df.join(count_df, on=["chrom", "pos"], how="left")
    
    return df


def count_snp_alleles(bam: AlignmentFile, chrom: str, snp_list: Iterable[Tuple[int, str, str]]) -> List[Tuple[str, int, int, int, int]]:
    """
    Compute allele counts for SNP positions in a given chromosome.

    This helper function iterates over an iterable of SNP tuples (each consisting of a position,
    reference allele, and alternate allele) for the specified chromosome. For each SNP, it fetches
    reads from the BAM file and counts the number of reads supporting the reference allele, alternate allele,
    or other bases. Each read is only counted once.

    Parameters
    ----------
    bam : pysam.libcalignmentfile.AlignmentFile
        An open BAM file object.
    chrom : str
        The chromosome name.
    snp_list : iterable of tuple
        An iterable of SNP tuples in the format (pos, ref, alt).

    Returns
    -------
    list of tuple
        A list of tuples in the format (chrom, pos, ref_count, alt_count, other_count).

    Notes
    -----
    The following commented code is preserved:
    
        # read_set = set()
        
        # TODO Update with binary search
        # Found no longer need to loop
    """
    read_set = set()
    allele_counts: List[Tuple[str, int, int, int, int]] = []

    for pos, ref, alt in snp_list:
        # read_set = set()
        ref_count, alt_count, other_count = 0, 0, 0

        # Got make sure read is not double counted
        for read in bam.fetch(chrom, pos - 1, pos):
            # If already counted allele
            if read.query_name in read_set:
                continue
            
            read_set.add(read.query_name)
            
            seq = read.query_sequence
            
            for qpos, refpos in read.get_aligned_pairs(True):
                # TODO Update with binary search
                if refpos == pos - 1:
                    if seq[qpos] == ref:
                        ref_count += 1
                    elif seq[qpos] == alt:
                        alt_count += 1
                    else:
                        other_count += 1
                    # Found no longer need to loop
                    break
        
        allele_counts.append((chrom, pos, ref_count, alt_count, other_count))
                
    return allele_counts
