"""
Author: Aaron Ho
Python Version: 3.8
"""

# Default Python package Imports
from pathlib import Path
from typing import List, Optional, Tuple

# External package imports
import pysam
import pandas as pd
from pysam import VariantFile
from pybedtools import BedTool


def write_sample_snp(in_file: str, in_sample: str, out_dir: str) -> None:
    """
    Filter heterozygous SNPs by sample and write to a new VCF file.

    This function reads a VCF file, subsets the samples to the specified sample,
    and filters for heterozygous SNPs (only considering SNPs with single nucleotide alleles).
    Filtered records are written to a new VCF file in the specified output directory.

    Parameters
    ----------
    in_file : str
        Path to the input VCF file.
    in_sample : str
        Name of the sample column in the VCF to check the genotype (GT).
    out_dir : str
        Output directory where the filtered VCF file will be written.

    Returns
    -------
    None

    Examples
    --------
    >>> write_sample_snp("input.vcf", "Sample1", "output_directory")
    """
    vcf = VariantFile(in_file)
    vcf.subset_samples([in_sample])
    
    out_vcf = VariantFile(str(Path(out_dir) / "filter.vcf"), "w", header=vcf.header)

    vcf_data = vcf.fetch()

    for record in vcf_data:
        if ((len(record.ref) == 1) and (len(record.alts) == 1) and (len(record.alts[0]) == 1)
                and (((record.samples[in_sample]['GT'][0] == 0) and (record.samples[in_sample]['GT'][1] == 1))
                     or ((record.samples[in_sample]['GT'][0] == 1) and (record.samples[in_sample]['GT'][1] == 0)))):
            out_vcf.write(record)

    print("Created Filtered VCF")


def write_filter_gtf(gtf_file: str, feature: Optional[List[str]], out_dir: str) -> None:
    """
    Filter a GTF file by feature and write the filtered data to a new file.

    This function reads a GTF file into a pandas DataFrame, filters the rows based on the
    provided feature(s) (if any), and writes the resulting DataFrame to a new GTF file in the
    specified output directory.

    Parameters
    ----------
    gtf_file : str
        Path to the input GTF file.
    feature : list or None
        Feature or list of features to filter by. If None, no filtering is performed.
    out_dir : str
        Output directory where the filtered GTF file will be written.

    Returns
    -------
    None

    Examples
    --------
    >>> write_filter_gtf("input.gtf", ["exon", "CDS"], "output_directory")
    """
    df = pd.read_csv(
        gtf_file, sep="\t", header=None,
        names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
        dtype=object
    )

    if feature is not None:
        df = df.loc[df["feature"].isin(feature)]

    if out_dir is not None:
        df.to_csv(str(Path(out_dir) / "filter.gtf"), sep="\t", header=False, index=False)
        print(f"GTF filtered by feature")


def intersect_snp(vcf_file: str, region_file: str, out_dir: str) -> None:
    """
    Retrieve SNPs that intersect given regions and write the output to a file.

    This function uses pybedtools to intersect a (filtered) VCF file with a region file (BED, Peaks, or GTF)
    and writes the resulting intersections to an output file in the specified directory.

    Parameters
    ----------
    vcf_file : str
        Path to the (filtered) VCF file.
    region_file : str
        Path to the region file (e.g., BED, Peaks, GTF) used for the intersection.
    out_dir : str
        Output directory where the intersected file will be written.

    Returns
    -------
    None

    Examples
    --------
    >>> intersect_snp("filter.vcf", "regions.bed", "output_directory")
    """
    a = BedTool(vcf_file)
    b = BedTool(region_file)

    a.intersect(b, wb=True, output=str(Path(out_dir) / "intersect.bed"))

    print("Created Intersection File")


def parse_intersect_df(intersect_file: str) -> pd.DataFrame:
    """
    Parse an intersection file and create a DataFrame of SNPs that intersect regions.

    This function reads the intersection file (created by intersect_snp) into a pandas DataFrame,
    renames the columns appropriately, and creates a unique identifier for each peak.

    Parameters
    ----------
    intersect_file : str
        Path to the intersection file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'chrom', 'pos', 'ref', 'alt', 'peak'.
        Each row corresponds to a SNP that intersects a region.

    Examples
    --------
    >>> df = parse_intersect_df("intersect.bed")
    """
    df = pd.read_csv(
        intersect_file, sep="\t", header=None, usecols=[0, 1, 3, 4, 10, 11, 12],
        dtype={11: str, 12: str}
    )
    df.columns = ["chrom", "pos", "ref", "alt", "peak_chrom", "peak_start", "peak_end"]
    df["peak"] = df["peak_chrom"] + "_" + df["peak_start"] + "_" + df["peak_end"]

    return_df = df[["chrom", "pos", "ref", "alt", "peak"]].drop_duplicates().reset_index(drop=True)

    print("SNP DF Created")
    return return_df


def parse_gene_df(intersect_file: str) -> pd.DataFrame:
    """
    Parse an intersection file and create a DataFrame of SNPs intersecting genes.

    This function reads the intersection file, extracts relevant columns, and attempts to
    extract gene names from the attributes column. It returns a DataFrame that includes gene
    names along with SNP information.

    Parameters
    ----------
    intersect_file : str
        Path to the intersection file created by intersect_snp.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'chrom', 'pos', 'ref', 'alt', 'feature', 'genes'.
        Each row corresponds to a SNP intersecting a gene region.

    Examples
    --------
    >>> df = parse_gene_df("intersect.bed")
    """
    df = pd.read_csv(
        intersect_file, sep="\t", header=None, usecols=[0, 1, 3, 4, 12, 18]
    )
    df.columns = ["chrom", "pos", "ref", "alt", "feature", "attributes"]

    df["genes"] = df["attributes"].str.extract(r'(?<=name\s)(.*?);')
    df["genes"] = df["genes"].str.strip('"')

    return_df = df[["chrom", "pos", "ref", "alt", "feature", "genes"]].drop_duplicates().reset_index(drop=True)

    print("SNP DF Created")
    return return_df


def process_bam(bam_file: str, region_file: str, out_dir: str) -> None:
    """
    Filter a BAM file to remove reads not overlapping regions of interest.

    This function uses pysam to filter a BAM file based on regions provided in a region file
    (BED, Peaks, or GTF). The filtered BAM is sorted and indexed, and the output is written
    to the specified output directory.

    Parameters
    ----------
    bam_file : str
        Path to the input BAM file.
    region_file : str
        Path to the region file (e.g., BED, Peaks, GTF) defining regions of interest.
    out_dir : str
        Output directory where the filtered BAM file and its index will be written.

    Returns
    -------
    None

    Examples
    --------
    >>> process_bam("input.bam", "regions.bed", "output_directory")
    """
    out_bam = Path(out_dir) / "filter.bam"
    sort_out = Path(out_dir) / "filter.sort.bam"

    print("Filtering reads that overlap regions of interest")
    pysam.view("-L", str(region_file), "-o", str(out_bam), str(bam_file), catch_stdout=False)
    pysam.sort(str(out_bam), "-o", str(sort_out), catch_stdout=False)
    pysam.index(str(sort_out), catch_stdout=False)

    print("Bam file filtered!")
