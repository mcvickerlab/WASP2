"""Data filtering utilities for allele-specific analysis.

Functions for filtering VCF, GTF, BAM files and creating intersection files.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import pysam
from pybedtools import BedTool
from pysam import VariantFile

if TYPE_CHECKING:
    from collections.abc import Sequence

logger = logging.getLogger(__name__)


def write_sample_snp(in_file: str | Path, in_sample: str, out_dir: str | Path) -> None:
    """
    Filters heterozygous SNP's by sample and writes to new VCF

    :param str in_file: Path to VCF file
    :param str in_sample: Name of sample column in VCF to check GT
    :param str out_dir: Name of output directory to write filtered VCF
    """
    vcf = VariantFile(str(in_file))
    vcf.subset_samples([in_sample])

    out_vcf = VariantFile(str(Path(out_dir) / "filter.vcf"), "w", header=vcf.header)

    vcf_data = vcf.fetch()

    for record in vcf_data:
        alts = record.alts
        ref = record.ref
        if (
            alts is not None
            and ref is not None
            and (len(ref) == 1)
            and (len(alts) == 1)
            and (len(alts[0]) == 1)
            and (
                (
                    (record.samples[in_sample]["GT"][0] == 0)
                    and (record.samples[in_sample]["GT"][1] == 1)
                )
                or (
                    (record.samples[in_sample]["GT"][0] == 1)
                    and (record.samples[in_sample]["GT"][1] == 0)
                )
            )
        ):
            out_vcf.write(record)

    logger.info("Created filtered VCF")


def write_filter_gtf(
    gtf_file: str | Path,
    feature: Sequence[str] | None,
    out_dir: str | Path | None,
) -> None:
    """Filter GTF file by feature type.

    Parameters
    ----------
    gtf_file : str | Path
        Path to GTF file.
    feature : Sequence[str] | None
        Feature types to keep (e.g., ['gene', 'exon']).
    out_dir : str | Path | None
        Output directory for filtered GTF.
    """
    df = pd.read_csv(
        gtf_file,
        sep="\t",
        header=None,
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
        dtype=object,
    )

    if feature is not None:
        df = df.loc[df["feature"].isin(feature)]

    if out_dir is not None:
        df.to_csv(str(Path(out_dir) / "filter.gtf"), sep="\t", header=False, index=False)
        logger.info("GTF filtered by feature")


def intersect_snp(vcf_file: str | Path, region_file: str | Path, out_dir: str | Path) -> None:
    """
    Retrieves SNP's that intersect regions

    :param str vcf_file: Path to (Filtered) VCF file
    :param str region_file: Path to region file (BED, Peaks, GTF)
    :param str out_dir: Name of output directory to write intersected VCF
    """
    a = BedTool(vcf_file)
    b = BedTool(region_file)

    a.intersect(b, wb=True, output=str(Path(out_dir) / "intersect.bed"))

    logger.info("Created intersection file")


def parse_intersect_df(intersect_file: str | Path) -> pd.DataFrame:
    """
    Parses intersection file and creates Dataframe

    :param intersect_file: Intersection file created by intersect_snp()
    :return DataFrame: Dataframe with SNP's that intersect regions
    """
    df = pd.read_csv(
        intersect_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 3, 4, 10, 11, 12],
        dtype={11: str, 12: str},
    )
    df.columns = ["chrom", "pos", "ref", "alt", "peak_chrom", "peak_start", "peak_end"]
    df["peak"] = df["peak_chrom"] + "_" + df["peak_start"] + "_" + df["peak_end"]

    return_df = df[["chrom", "pos", "ref", "alt", "peak"]].drop_duplicates().reset_index(drop=True)

    logger.info("SNP DataFrame created")
    return return_df


def parse_gene_df(intersect_file: str | Path) -> pd.DataFrame:
    """
    Parses intersection file and creates Dataframe
    Returns gene names

    :param intersect_file: Intersection file created by intersect_snp()
    :return DataFrame: Dataframe with SNP's that intersect regions
    """
    df = pd.read_csv(intersect_file, sep="\t", header=None, usecols=[0, 1, 3, 4, 12, 18])
    df.columns = ["chrom", "pos", "ref", "alt", "feature", "attributes"]

    df["genes"] = df["attributes"].str.extract(r"(?<=name\s)(.*?);")
    df["genes"] = df["genes"].str.strip('"')

    return_df = (
        df[["chrom", "pos", "ref", "alt", "feature", "genes"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    logger.info("SNP DataFrame created")
    return return_df


def process_bam(bam_file: str | Path, region_file: str | Path, out_dir: str | Path) -> None:
    """
    Filter bam file to remove reads not overlapping regions of interest

    :param str bam_file: Path to BAM file
    :param str region_file: Path to region file (BED, Peaks, GTF)
    :param str out_dir: Path to output directory of filtered BAM
    """
    out_bam = Path(out_dir) / "filter.bam"
    sort_out = Path(out_dir) / "filter.sort.bam"

    logger.info("Filtering reads that overlap regions of interest")
    pysam.view("-L", str(region_file), "-o", str(out_bam), str(bam_file), catch_stdout=False)
    pysam.sort(str(out_bam), "-o", str(sort_out), catch_stdout=False)
    pysam.index(str(sort_out), catch_stdout=False)

    logger.info("BAM file filtered")
