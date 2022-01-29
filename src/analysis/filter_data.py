"""
Author: Aaron Ho
Python Version: 3.8
"""

# Default Python package Imports
from pathlib import Path

# External package imports
import pysam
import pandas as pd
from pysam import VariantFile
from pybedtools import BedTool


def write_sample_snp(in_file, in_sample, out_dir):
    """
    Filters heterozygous SNP's by sample and writes to new VCF

    :param str in_file: Path to VCF file
    :param str in_sample: Name of sample column in VCF to check GT
    :param str out_dir: Name of output directory to write filtered VCF
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


def write_filter_gtf(gtf_file, feature, out_dir):
    df = pd.read_csv(gtf_file, sep="\t", header=None,
     names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
      dtype=object)

    if feature is not None:
        df = df.loc[df["feature"].isin(feature)]

    if out_dir is not None:
        df.to_csv(str(Path(out_dir) / "filter.gtf"), sep="\t", header=False, index=False)
        print(f"GTF filtered by feature")


def intersect_snp(vcf_file, region_file, out_dir):
    """
    Retrieves SNP's that intersect regions

    :param str vcf_file: Path to (Filtered) VCF file
    :param str region_file: Path to region file (BED, Peaks, GTF)
    :param str out_dir: Name of output directory to write intersected VCF
    """
    a = BedTool(vcf_file)
    b = BedTool(region_file)

    a.intersect(b, wb=True, output=str(Path(out_dir) / "intersect.bed"))

    print("Created Intersection File")


def parse_intersect_df(intersect_file):
    """
    Parses intersection file and creates Dataframe

    :param intersect_file: Intersection file created by intersect_snp()
    :return DataFrame: Dataframe with SNP's that intersect regions
    """
    df = pd.read_csv(intersect_file, sep="\t", header=None, usecols=[0, 1, 3, 4, 10, 11, 12], dtype={11: str, 12: str})
    df.columns = ["chrom", "pos", "ref", "alt", "peak_chrom", "peak_start", "peak_end"]
    df["peak"] = df["peak_chrom"] + "_" + df["peak_start"] + "_" + df["peak_end"]

    return_df = df[["chrom", "pos", "ref", "alt", "peak"]].drop_duplicates().reset_index(drop=True)

    print("SNP DF Created")
    return return_df


def parse_gene_df(intersect_file):
    """
    Parses intersection file and creates Dataframe
    Returns gene names

    :param intersect_file: Intersection file created by intersect_snp()
    :return DataFrame: Dataframe with SNP's that intersect regions
    """
    df = pd.read_csv(intersect_file, sep="\t", header=None, usecols=[0, 1, 3, 4, 12, 18])
    df.columns = ["chrom", "pos", "ref", "alt", "feature", "attributes"]

    df["genes"] = df["attributes"].str.extract(r'(?<=name\s)(.*?);')
    df["genes"] = df["genes"].str.strip('"')

    return_df = df[["chrom", "pos", "ref", "alt", "feature", "genes"]].drop_duplicates().reset_index(drop=True)

    print("SNP DF Created")
    return return_df


def process_bam(bam_file, region_file, out_dir):
    """
    Filter bam file to remove reads not overlapping regions of interest

    :param str bam_file: Path to BAM file
    :param str region_file: Path to region file (BED, Peaks, GTF)
    :param str out_dir: Path to output directory of filtered BAM
    """

    out_bam = Path(out_dir) / "filter.bam"
    sort_out = Path(out_dir) / "filter.sort.bam"

    print("Filtering reads that overlap regions of interest")
    pysam.view("-L", str(region_file), "-o", str(out_bam), str(bam_file), catch_stdout=False)
    pysam.sort(str(out_bam), "-o", str(sort_out), catch_stdout=False)
    pysam.index(str(sort_out), catch_stdout=False)

    print("Bam file filtered!")
