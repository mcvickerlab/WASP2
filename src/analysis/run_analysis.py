"""
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
from pathlib import Path
import argparse
import time
import tempfile

# External package imports
import pandas as pd

# Local script imports
from filter_data import write_sample_snp, intersect_snp, parse_intersect_df, process_bam
from count_alleles import make_count_df
from count_alleles_sc import make_count_df_sc
from as_analysis import get_imbalance, get_imbalance_sc


def preprocess_data(in_bam, in_vcf, in_region, in_sample, nofilt, out_dir):

    if not nofilt:
        process_bam(in_bam, in_region, out_dir)

    write_sample_snp(in_vcf, in_sample, out_dir)

    intersect_snp(str(Path(out_dir) / "filter.vcf"), in_region, out_dir)

    intersect_df = parse_intersect_df(str(Path(out_dir) / "intersect.bed"))

    return intersect_df
    

def parse_counting(in_bam, in_vcf, in_region, in_sample, out_dir, nofilt=False, temp_loc=None):

    start = time.time()
    if temp_loc is None:
        with tempfile.TemporaryDirectory() as tmpdir:

            intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, nofilt, tmpdir)

            if nofilt is True:
                df = make_count_df(in_bam, intersect_df)
            else:
                df = make_count_df(str(Path(tmpdir) / "filter.sort.bam"), intersect_df)

            df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)

    else:
        intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, nofilt, temp_loc)

        if nofilt is True:
            df = make_count_df(in_bam, intersect_df)
        else:
            df = make_count_df(str(Path(temp_loc) / "filter.sort.bam"), intersect_df)
        
        df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)
    
    end = time.time()
    print(f"Counting pipeline completed in {end - start} seconds!")


def parse_counting_sc(in_bam, in_vcf, in_region, in_sample, in_barcodes, out_dir, nofilt=False, temp_loc=None):
    bc_series = pd.read_csv(in_barcodes, sep="\t", header=None, index_col=0, names=["barcodes", "group"], squeeze=True)

    start = time.time()
    if temp_loc is None:
        with tempfile.TemporaryDirectory() as tmpdir:

            intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, nofilt, tmpdir)

            if nofilt is True:
                df = make_count_df_sc(in_bam, intersect_df, bc_series)
            else:
                df = make_count_df_sc(str(Path(tmpdir) / "filter.sort.bam"), intersect_df, bc_series)

            df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)

    else:
        intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, nofilt, temp_loc)

        if nofilt is True:
            df = make_count_df_sc(in_bam, intersect_df, bc_series)
        else:
            df = make_count_df_sc(str(Path(temp_loc) / "filter.sort.bam"), intersect_df, bc_series)
        
        df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)
    
    end = time.time()
    print(f"Counting pipeline completed in {end - start} seconds!")


def main():
    parent_parser = argparse.ArgumentParser(add_help=False)
    seq_type = parent_parser.add_mutually_exclusive_group() # Denote rna-seq vs atac-seq
    seq_type.add_argument("--rna", action='store_true', help="TODO") # TODO make RNA-SEQ work with software
    seq_type.add_argument("--atac", action='store_true', help="default")

    parent_parser.add_argument("-sc", "--singlecell", action="store_true", help="Single Cell Option")


    parser = argparse.ArgumentParser()

    # Subparser options
    subparser = parser.add_subparsers(dest="command")

    count_parser = subparser.add_parser("count", parents=[parent_parser])
    count_parser.add_argument("-a", "--alignment", required=True, type=str,help="Alignment BAM File")
    count_parser.add_argument("-g", "--genotypes", required=True, type=str, help="Genotype VCF File")
    count_parser.add_argument("-s", "--sample", required=True, type=str, help="Sample name in VCF")
    count_parser.add_argument("-r", "--regions", required=True, type=str, help="Genes or Peak File in BED, narrowPeak, or GTF Format")
    count_parser.add_argument("-b", "--barcodes", type=str, help="Two row TSV mapping barcode to cell-group/cluster")
    count_parser.add_argument("-o", "--output", type=str, help="Output Directory", default=str(Path.cwd()))

    count_parser.add_argument("--nofilt", action='store_true') # TODO filtering options
    count_parser.add_argument("--keeptemps", type=str, nargs="?", const=str(Path.cwd())) # TODO allow for intermediate files to be saved

   

    analysis_parser = subparser.add_parser("analysis", parents=[parent_parser])
    analysis_parser.add_argument("counts")
    analysis_parser.add_argument("--min", type=int, help="Minimum allele count for analysis", default=10)
    analysis_parser.add_argument("-m", "--model", type=str, choices=["single", "linear", "binomial"], help="Analysis Model", default="single")
    analysis_parser.add_argument("-o", "--output", type=str, help="Output Directory", default=str(Path.cwd()))


    args = parser.parse_args()

    # print(args)

    if args.singlecell: # Single Cell Data
        print("Single Cell Analysis")
        if (args.command == "count") and not (args.barcodes):
            parser.error("count --singlecell requires --barcodes")

        elif (args.command == "count") and args.barcodes:
            parse_counting_sc(args.alignment, args.genotypes, args.regions, args.sample, args.barcodes, args.output, nofilt=args.nofilt, temp_loc=args.keeptemps)
        
        elif args.command == "analysis":
            get_imbalance_sc(args.counts, args.min, args.model, args.output)


    else: # Bulk processing
        print("Bulk Analysis")
        if (args.command == "count") and not (args.barcodes):
            parse_counting(args.alignment, args.genotypes, args.regions, args.sample, args.output, nofilt=args.nofilt, temp_loc=args.keeptemps)

        elif (args.command == "count") and args.barcodes:
            parser.error("Single-Cell barcodes given for bulk analysis. Remove --barcode input, or instead run with --singlecell option")

        elif args.command == "analysis":
            get_imbalance(args.counts, args.min, args.model, args.output)


if __name__ == '__main__':
    main()
