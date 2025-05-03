from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Tuple

import timeit
import numpy as np
import pandas as pd
import polars as pl
import anndata as ad

from scipy.sparse import csr_matrix
from pysam.libcalignmentfile import AlignmentFile

# Local imports
from wasp2.counting.count_alleles import find_read_aln_pos


class CountStatsSC:
    """
    Class to hold mutable and persistent single-cell count statistics.

    This class uses defaultdicts to maintain counts and metadata across chromosomes,
    including reference, alternate, and other allele counts, as well as various read
    and SNP statistics.
    """

    def __init__(self) -> None:
        """
        Initialize all count and metadata attributes as defaultdicts of int.
        """
        self.ref_count: defaultdict = defaultdict(int)
        self.alt_count: defaultdict = defaultdict(int)
        self.other_count: defaultdict = defaultdict(int)
        
        # Keep track of metadata
        
        # Number
        self.num_snps: defaultdict = defaultdict(int)
        self.num_barcodes: defaultdict = defaultdict(int)
        self.reads_counted: defaultdict = defaultdict(int)
        
        # Reads that were not counted
        self.reads_skipped_no_barcode: defaultdict = defaultdict(int)
        self.reads_skipped_barcode_no_index: defaultdict = defaultdict(int)
        self.reads_skipped_prev_counted: defaultdict = defaultdict(int)

    def stats_to_df(self) -> pd.DataFrame:
        """
        Convert count statistics to a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'chrom', 'num_snps', 'num_barcodes', 'reads_counted',
            'reads_skipped_no_barcode', 'reads_skipped_barcode_no_index', and 'reads_skipped_prev_counted'.
        """
        stat_attributes = [
            "num_snps", "num_barcodes", "reads_counted",
            "reads_skipped_no_barcode",
            "reads_skipped_barcode_no_index",
            "reads_skipped_prev_counted"
        ]
        
        stat_df = pd.DataFrame(
            {key: getattr(self, key) for key in stat_attributes}
        ).reset_index(names="chrom")
        
        return stat_df


def make_count_matrix(
    bam_file: str,
    df: pl.DataFrame,
    bc_dict: Dict[Any, int],
    include_samples: Optional[List[str]] = None,
    include_features: Optional[List[str]] = None,
) -> ad.AnnData:
    """
    Create a sparse count matrix and annotate an AnnData object with allele counts.

    This function extracts SNP information from a Polars DataFrame, counts alleles using a BAM file,
    and creates sparse matrices for reference, alternate, and other allele counts. These matrices are then
    combined into an AnnData object which is further annotated with metadata and, optionally, sample and
    feature information.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    df : polars.DataFrame
        DataFrame containing intersection data with a "chrom" column.
    bc_dict : dict
        Dictionary mapping barcodes to column indices.
    include_samples : list, optional
        List of sample names to include in the annotations.
    include_features : list, optional
        List of feature names to include (e.g., 'region' or gene information).

    Returns
    -------
    anndata.AnnData
        Annotated AnnData object with allele count layers and metadata.
    
    Notes
    -----
    The following commented code is preserved:
    
        # chrom_list = chrom_list[:3] # Testing purposes
        # Add genotypes annotations
        # Maybe do this automatically and parse feature col instead?
        # Might be more memory efficient to use pandas index instead...
        # region_snp_dict = dict(
        #     df.join(snp_df, on=["chrom", "pos", "ref", "alt"], how="left"
        #            ).group_by("region", maintain_order=True
        #                      ).agg("index").iter_rows()
        # )
        # adata.uns["region_snps"] = region_snp_dict
    """
    chrom_list: List[str] = df.get_column("chrom").unique(maintain_order=True)
    # chrom_list = chrom_list[:3] # Testing purposes

    # Add genotypes annotations
    # Maybe do this automatically and parse feature col instead?
    snp_df_cols: List[str] = ["chrom", "pos", "ref", "alt"]
    if include_samples is not None:
        snp_df_cols.extend(include_samples)
    
    # Might be more memory efficient to use pandas index instead...
    snp_df = df.select(snp_df_cols).unique(maintain_order=True).with_row_index()
    
    sc_counts = CountStatsSC()  # Class that holds total count data

    with AlignmentFile(bam_file, "rb") as bam:
        for chrom in chrom_list:
            chrom_df = snp_df.filter(pl.col("chrom") == chrom)
            start = timeit.default_timer()
            
            try:
                count_bc_snp_alleles(
                    bam=bam,
                    bc_dict=bc_dict,
                    chrom=chrom,
                    snp_list=chrom_df.select(["index", "pos", "ref", "alt"]).iter_rows(),
                    sc_counts=sc_counts,
                )
            except ValueError:
                print(f"Skipping {chrom}: Contig not found!")
            else:
                print(f"{chrom}: Counted {chrom_df.height} SNP's in {timeit.default_timer() - start:.2f} seconds!")

    # Create sparse matrices
    # sparse array is recommended...but doesnt work with adata
    sparse_ref = csr_matrix(
        (list(sc_counts.ref_count.values()),
         list(zip(*sc_counts.ref_count.keys()))),
        shape=(snp_df.shape[0], len(bc_dict)),
        dtype=np.uint8,
    )
    
    sparse_alt = csr_matrix(
        (list(sc_counts.alt_count.values()),
         list(zip(*sc_counts.alt_count.keys()))),
        shape=(snp_df.shape[0], len(bc_dict)),
        dtype=np.uint8,
    )
    
    sparse_other = csr_matrix(
        (list(sc_counts.other_count.values()),
         list(zip(*sc_counts.other_count.keys()))),
        shape=(snp_df.shape[0], len(bc_dict)),
        dtype=np.uint8,
    )
    
    # Create anndata With total as X
    adata = ad.AnnData(
        X=sparse_ref + sparse_alt + sparse_other,
        layers={
            "ref": sparse_ref,
            "alt": sparse_alt,
            "other": sparse_other,
        },
    )

    # Annotate adata: Figure out what to add to adata here vs later
    adata.obs = snp_df.to_pandas()  # Maybe just switch to pandas? Should i set no copy?
    adata.obs["ref_count"] = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1
    adata.obs["alt_count"] = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1
    
    # Add barcode names
    adata.var_names = list(bc_dict.keys())
    
    # Add genotypes to anndata
    if include_samples is not None:
        adata.uns["samples"] = include_samples

    # TODO: Allow for other features besides 'region' using include_features
    # Could be case of no features, or feature is gene
    if "region" in df.columns:
        # Get unique snps and associated regions
        # Create dict during analysis step instead
        adata.uns["feature"] = df.join(snp_df,
                                       on=["chrom", "pos", "ref", "alt"],
                                       how="left").select(
                                           ["region", "index"]).to_pandas()
        
        # region_snp_dict = dict(
        #     df.join(snp_df, on=["chrom", "pos", "ref", "alt"], how="left"
        #            ).group_by("region", maintain_order=True
        #                      ).agg("index").iter_rows()
        # )
        # adata.uns["region_snps"] = region_snp_dict

    # Write out count stats
    adata.uns["count_stats"] = sc_counts.stats_to_df()
    
    return adata


def count_bc_snp_alleles(
    bam: AlignmentFile,
    bc_dict: Dict[Any, int],
    chrom: str,
    snp_list: Iterable[Tuple[Any, Any, Any, Any]],
    sc_counts: CountStatsSC,
) -> None:
    """
    Count SNP alleles for a given chromosome and update single-cell count statistics.

    This function fetches reads from a BAM file for a specific chromosome, iterates over SNP entries,
    and updates a CountStatsSC instance with the counts for reference, alternate, and other alleles.
    It uses barcode mapping (bc_dict) to assign counts to the appropriate column index and tracks various
    statistics such as the number of reads counted and SNPs processed.

    Parameters
    ----------
    bam : pysam.libcalignmentfile.AlignmentFile
        An open BAM file.
    bc_dict : dict
        Dictionary mapping barcodes to column indices.
    chrom : str
        Chromosome name.
    snp_list : iterable of tuple
        An iterable of SNP tuples in the format (index, pos, ref, alt).
    sc_counts : CountStatsSC
        An instance of CountStatsSC to update count statistics.

    Returns
    -------
    None
    """
    read_set = set()  # Keep track of reads seen
    bc_set = set()

    for idx, pos, ref, alt in snp_list:
        for read in bam.fetch(chrom, pos - 1, pos):
            # If already counted allele or pair in read
            if read.query_name in read_set:
                sc_counts.reads_skipped_prev_counted[chrom] += 1
                continue

            # Check if there is a read barcode
            try:
                read_bc = read.get_tag("CB")
            except KeyError:
                sc_counts.reads_skipped_no_barcode[chrom] += 1
                continue

            # If barcode not in dict
            if read_bc not in bc_dict:
                sc_counts.reads_skipped_barcode_no_index[chrom] += 1
                continue

            seq = read.query_sequence

            # TEST binary search
            qpos = find_read_aln_pos(read, pos - 1)

            try:
                if seq[qpos] == ref:
                    sc_counts.ref_count[(idx, bc_dict[read_bc])] += 1
                elif seq[qpos] == alt:
                    sc_counts.alt_count[(idx, bc_dict[read_bc])] += 1
                else:
                    sc_counts.other_count[(idx, bc_dict[read_bc])] += 1
            except TypeError:
                continue
            else:
                read_set.add(read.query_name)
                bc_set.add(read_bc)
                sc_counts.reads_counted[chrom] += 1

        sc_counts.num_snps[chrom] += 1  # Put here in case of error
    sc_counts.num_barcodes[chrom] = len(bc_set)  # Add unique barcodes
