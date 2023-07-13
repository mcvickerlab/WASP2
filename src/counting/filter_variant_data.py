import sys
import timeit
import subprocess
import warnings


from pathlib import Path

import numpy as np
import polars as pl

# same as in mapping...should create unified utils
def vcf_to_bed(vcf_file, out_bed, samples=None, drop_gt=False):
    
    # Base commands
    view_cmd = ["bcftools", "view", str(vcf_file),
                "-m2", "-M2", "-v", "snps", "-Ou"
               ]

    query_cmd = ["bcftools", "query",
                 "-o", str(out_bed),
                 "-f"]
    
    # Parse based on num samps
    if samples is None:
        
        # 0 samps, no GTs
        view_cmd.append("--drop-genotypes")
        query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")
    else:
        
        # Samples
        samples_arg = ",".join(samples)
        num_samples = len(samples)
        
        if num_samples > 1:
            # Multisamp
            view_cmd.extend(["-s", samples_arg,
                             "--min-ac", "1",
                             "--max-ac", str((num_samples * 2) - 1)])
        else:
            # Single Samp
            view_cmd.extend(["-s", samples_arg, "--genotype", "het"])
        
        # If we include GT
        if drop_gt:
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT\n")
        else:
            query_cmd.append("%CHROM\t%POS0\t%END\t%REF\t%ALT[\t%TGT]\n")
    
    # Run Subprocess
    view_process = subprocess.run(view_cmd, stdout=subprocess.PIPE, check=True)
    query_process = subprocess.run(query_cmd, input=view_process.stdout, check=True)
    
    # return out_bed


# Perform intersection
def intersect_vcf_region(vcf_file, region_file, out_file):
    
    # Parse region file before or after???
    intersect_cmd = ["bedtools", "intersect" , "-a", str(vcf_file),
                     "-b", str(region_file), "-wb"]
    

    # write intersect out
    with open(out_file, "w") as file:
        intersect_process = subprocess.run(intersect_cmd, stdout=file, check=True)


# Convert Intersect file to df
def parse_intersect_region(intersect_file, use_region_names=False):
    df = pl.scan_csv(intersect_file, separator="\t", has_header=False)
    
    # If we need to use coords as name
    use_coords = False
    
    # No regions, only variants
    if len(df.columns) <= 5:
        # No regions, only variants
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4]]
        new_cols = ["chrom", "pos", "ref", "alt"]
    
    elif use_region_names and len(df.columns) >= 9:
        # Use included names in region file
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, 8]]
        new_cols = ["chrom", "pos", "ref", "alt", "region"]
        
    elif len(df.columns) >= 8:
        # Either no names included or use coords instead
        subset_cols = [df.columns[i] for i in [0, 2, 3, 4, 5, 6, 7]]
        new_cols = ["chrom", "pos", "ref", "alt",
                    "region_chrom", "region_start", "region_end"]
        use_coords = True

    else:
        # CHANGE TO RAISE ERROR
        print("COULD NOT RECOGNIZE FORMAT OR WRONG NUMBER OF COLS")
        return
    
    
    # Parse dataframe columns
    rename_cols = {old_col: new_col for old_col, new_col in zip(subset_cols, new_cols)}
    df = df.select(subset_cols).rename(
        rename_cols).with_columns(pl.col("pos").cast(pl.UInt32))
    
    # Create coords
    if use_coords:
        df = df.with_columns(
            pl.concat_str(
                [pl.col(i) for i in new_cols[-3::]],
                separator="_"
            ).alias("region")
        ).select(
            ["chrom", "pos", "ref", "alt", "region"])
    
    return df.unique(maintain_order=True).collect()

