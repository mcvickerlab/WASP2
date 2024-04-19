import timeit
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import polars as pl
import anndata as ad

from scipy.sparse import csr_matrix
from pysam.libcalignmentfile import AlignmentFile

# Local imports
from count_alleles import find_read_aln_pos


# Create class that holds mutable and persistent stats
class CountStatsSC:

    def __init__(self):
        
        self.ref_count = defaultdict(int)
        self.alt_count = defaultdict(int)
        self.other_count = defaultdict(int)
        
        # Keep track of metadata
        
        # Number 
        self.num_snps = defaultdict(int)
        self.num_barcodes = defaultdict(int)
        self.reads_counted = defaultdict(int)
        
        # Reads that were not counted
        self.reads_skipped_no_barcode = defaultdict(int)
        self.reads_skipped_barcode_no_index = defaultdict(int)
        self.reads_skipped_prev_counted = defaultdict(int)
    
    def stats_to_df(self):
        
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


# Create sparse count matrix
def make_count_matrix(bam_file, df, bc_dict,
                      include_samples=None,
                      include_features=None
                     ):
    
    chrom_list = df.get_column("chrom").unique(maintain_order=True)
    # chrom_list = chrom_list[:3] # Testing purposes


    # Add genotypes annotations
    # Maybe do this automatically and parse feature col instead?
    snp_df_cols = ["chrom", "pos", "ref", "alt"]
    if include_samples is not None:
        snp_df_cols.extend(include_samples)
    
    
    # Might be more memory efficient to use pandas index instead...
    snp_df = df.select(snp_df_cols).unique(maintain_order=True).with_row_index()
    
    sc_counts = CountStatsSC() # Class that holds total count data

    with AlignmentFile(bam_file, "rb") as bam:

        for chrom in chrom_list:

            chrom_df = snp_df.filter(pl.col("chrom") == chrom)

            start = timeit.default_timer()
            
            try:
                
                count_bc_snp_alleles(
                    bam=bam,
                    bc_dict=bc_dict,
                    chrom=chrom,
                    snp_list=chrom_df.select(
                        ["index", "pos", "ref", "alt"]).iter_rows(),
                    sc_counts=sc_counts
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
        dtype=np.uint8
    )
    
    sparse_alt = csr_matrix(
        (list(sc_counts.alt_count.values()),
         list(zip(*sc_counts.alt_count.keys()))),
        shape=(snp_df.shape[0], len(bc_dict)),
        dtype=np.uint8
    )
    
    sparse_other = csr_matrix(
        (list(sc_counts.other_count.values()),
         list(zip(*sc_counts.other_count.keys()))),
        shape=(snp_df.shape[0], len(bc_dict)),
        dtype=np.uint8
    )
    
    
    # Create anndata With total as X
    adata = ad.AnnData(
        X=sparse_ref+sparse_alt+sparse_other,
        layers={
            "ref": sparse_ref,
            "alt": sparse_alt,
            "other": sparse_other
        }
    )

    
    # Annotate adata: Figure out what to add to adata here vs later
    adata.obs = snp_df.to_pandas() # Maybe just switch to pandas? Should i set no copy?
    adata.obs["ref_count"] = adata.layers["ref"].sum(axis=1, dtype=np.uint16).T.A1
    adata.obs["alt_count"] = adata.layers["alt"].sum(axis=1, dtype=np.uint16).T.A1
    
    
    # Add barcode names
    adata.var_names = bc_dict.keys()
    
    # Add genotypes to anndata
    if include_samples is not None:
        adata.uns["samples"] = include_samples

    
    # TODO: Allow for other features besides 'region' using include_features
    # Could be case of no features, or feature is gene
    if "region" in df.columns:
        
        # Get unique snps and associated regions
        region_snp_dict = dict(
            df.join(snp_df, on=["chrom", "pos", "ref", "alt"], how="left"
                   ).group_by("region", maintain_order=True
                             ).agg("index").iter_rows()
        )

        adata.uns["region_snps"] = region_snp_dict

        
    # Write out count stats
    adata.uns["count_stats"] = sc_counts.stats_to_df()
    

    return adata

    
def count_bc_snp_alleles(bam, bc_dict, chrom, snp_list, sc_counts):
    
    read_set = set() # Keep track of reads seen
    bc_set = set()

    
    for idx, pos, ref, alt in snp_list:
        
        for read in bam.fetch(chrom, pos-1, pos):

            # If already counted allele or pair in read
            if read.query_name in read_set:
                sc_counts.reads_skipped_prev_counted[chrom]+=1
                continue

            # Check if there is a read barcode
            try:
                read_bc = read.get_tag("CB")
            except KeyError:
                sc_counts.reads_skipped_no_barcode[chrom]+=1
                continue

            # If barcode not in dict
            if read_bc not in bc_dict:
                sc_counts.reads_skipped_barcode_no_index[chrom]+=1
                continue

            seq = read.query_sequence

            # TEST binary search
            qpos = find_read_aln_pos(read, pos-1)

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
                sc_counts.reads_counted[chrom]+=1
        
        sc_counts.num_snps[chrom]+=1 # Put here in case of error
    sc_counts.num_barcodes[chrom]=len(bc_set) # Add unique barcodes
