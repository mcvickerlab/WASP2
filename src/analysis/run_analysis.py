"""
Author: Aaron Ho
Python Version: 3.9
"""

# Default Python package Imports
from pathlib import Path
from csv import DictReader
import argparse
import time
import tempfile
import re

# External package imports
import pandas as pd
from pandas.io import parsers

# Local script imports
from filter_data import write_sample_snp, intersect_snp, parse_intersect_df, parse_gene_df, process_bam, write_filter_gtf
from count_alleles import make_count_df
from count_alleles_sc import make_count_df_sc
from as_analysis import get_imbalance, get_imbalance_sc


def preprocess_data(in_bam, in_vcf, in_region, in_sample, stype, nofilt, out_dir, features=None):

    if (stype == "rna") and features is not None:
        write_filter_gtf(in_region, features, out_dir)
        in_region = str(Path(out_dir) / "filter.gtf")

    if not nofilt:
        process_bam(in_bam, in_region, out_dir)

    write_sample_snp(in_vcf, in_sample, out_dir)

    intersect_snp(str(Path(out_dir) / "filter.vcf"), in_region, out_dir)

    if stype == "rna":
        intersect_df = parse_gene_df(str(Path(out_dir) / "intersect.bed"))
    else:
        intersect_df = parse_intersect_df(str(Path(out_dir) / "intersect.bed"))


    return intersect_df
    

def parse_counting(in_bam, in_vcf, in_region, in_sample, out_dir, stype, nofilt=False, temp_loc=None, features=None):

    start = time.time()
    if temp_loc is None:
        with tempfile.TemporaryDirectory() as tmpdir:

            intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, stype, nofilt, tmpdir, features)

            if nofilt is True:
                df = make_count_df(in_bam, intersect_df)
            else:
                df = make_count_df(str(Path(tmpdir) / "filter.sort.bam"), intersect_df)

            df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)

    else:
        intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, stype, nofilt, temp_loc, features)

        if nofilt is True:
            df = make_count_df(in_bam, intersect_df)
        else:
            df = make_count_df(str(Path(temp_loc) / "filter.sort.bam"), intersect_df)
        
        df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)
    
    end = time.time()
    print(f"Counting pipeline completed in {end - start} seconds!")


def parse_counting_sc(in_bam, in_vcf, in_region, in_sample, in_barcodes, out_dir, stype, nofilt=False, temp_loc=None, features=None):
    bc_series = pd.read_csv(in_barcodes, sep="\t", header=None, index_col=0, names=["barcodes", "group"], squeeze=True)

    start = time.time()
    if temp_loc is None:
        with tempfile.TemporaryDirectory() as tmpdir:

            intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, stype, nofilt, tmpdir, features)

            if nofilt is True:
                df = make_count_df_sc(in_bam, intersect_df, bc_series)
            else:
                df = make_count_df_sc(str(Path(tmpdir) / "filter.sort.bam"), intersect_df, bc_series)

            df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)

    else:
        intersect_df = preprocess_data(in_bam, in_vcf, in_region, in_sample, stype, nofilt, temp_loc, features)

        if nofilt is True:
            df = make_count_df_sc(in_bam, intersect_df, bc_series)
        else:
            df = make_count_df_sc(str(Path(temp_loc) / "filter.sort.bam"), intersect_df, bc_series)
        
        df.to_csv(str(Path(out_dir) / "as_counts.tsv"), sep="\t", header=True, index=False)
    
    end = time.time()
    print(f"Counting pipeline completed in {end - start} seconds!")


def parse_analysis(count_file, min_count, model, out_dir, stype, features=None):

    if stype == "rna":
        df = pd.read_csv(count_file, sep="\t")

        if features:
            df = df.loc[df["feature"].isin(features)]
        
        feature_list = df["feature"].unique()

        for feat in feature_list:
            feat_df = df.loc[df["feature"] == feat]
            feat_df = feat_df.drop(columns=["feature"])

            get_imbalance(feat_df, min_count, model, out_dir, is_gene=True, feature=feat)

    else:
        get_imbalance(count_file, min_count, model, out_dir, is_gene=False)


def parse_analysis_sc(count_file, min_count, model, out_dir, stype, features=None):

    if stype == "rna":
        df = pd.read_csv(count_file, sep="\t")

        if features:
            df = df.loc[df["feature"].isin(features)]
        
        feature_list = df["feature"].unique()

        for feat in feature_list:
            feat_df = df.loc[df["feature"] == feat]
            feat_df = feat_df.drop(columns=["feature"])

            get_imbalance_sc(feat_df, min_count, model, out_dir, is_gene=True, feature=feat)
        
    else:
        get_imbalance_sc(count_file, min_count, model, out_dir, is_gene=False)


def validate_args(args):    # TODO Better parsing of valid files and inputs
    singlecell = args.singlecell
    command = args.command

    if command == "count":

        # If no directory given for --keeptemps, set to outdir
        if args.keeptemps == 0:
            args.keeptemps = args.output
        

        # Manual parsing of required inputs for single-cell
        if singlecell and not args.barcodes:
            return "count --singlecell requires --barcodes"
        elif not singlecell and args.barcodes:
            return "Single-Cell barcodes given for bulk analysis. Remove --barcode input, or instead run with --singlecell option"
        

        # Validate if rna-seq or ATAC-seq
        f_ext = "".join(Path(args.regions).suffixes)

        if re.search(r'\.(.*Peak|bed)(\.gz)?$', f_ext, re.I):
            file_stype = "atac"
        elif re.search(r'\.(gtf)(\.gz)?$', f_ext, re.I):
            file_stype = "rna"
        else:
            return "Unrecognized file-type for region file"

        if args.stype is None:
            args.stype = file_stype
        elif args.stype != file_stype:
            return f"--{args.stype} option selected, but {f_ext} file given"
        

        # Validate feature options
        if (args.stype != "rna") and (args.features is not None):
            return f"--features only valid using rna-seq data, and gtf annotations"

        elif (args.stype == "rna") and (args.features is None):
            args.features = ["transcript"]
        
        elif (args.stype == "rna") and not args.features:
            args.features = None
    
    elif command =="analysis":

        # Confirm input type
        if args.stype is None: # Use columns to determine rna or atac
            with open(args.counts, "r") as file:
                header=DictReader(file, delimiter="\t").fieldnames
            
            if "genes" in header:
                args.stype = "rna"
            else:
                args.stype = "atac"

        if (args.stype != "rna") and (args.features is not None):
            return f"--features only valid using rna-seq data"

    # TODO more sanity checks, as well as parsing different file types
    return args


def parse_cmd():
    parent_parser = argparse.ArgumentParser(add_help=False)
    seq_type = parent_parser.add_mutually_exclusive_group()  # Denote rna-seq vs atac-seq
    seq_type.add_argument("--rna", action='store_const', const="rna", dest="stype",
                          help="Explicitly denote RNA-seq data (Otherwise infers from input)")
    seq_type.add_argument("--atac", action='store_const', const="atac", dest="stype",
                          help="Explicitly denote ATAC-seq data (Otherwise infers from input)")

    parent_parser.add_argument("-sc", "--singlecell", action="store_true", help="Single Cell Option")

    parent_parser.add_argument("-ft", "--features", nargs="*", help="GTF features to analyze")

    parser = argparse.ArgumentParser()

    # Subparser options
    subparser = parser.add_subparsers(dest="command")

    count_parser = subparser.add_parser("count", parents=[parent_parser])
    count_parser.add_argument("-a", "--alignment", required=True, type=str, help="Alignment BAM File")
    count_parser.add_argument("-g", "--genotypes", required=True, type=str, help="Genotype VCF File")
    count_parser.add_argument("-s", "--sample", required=True, type=str, help="Sample name in VCF")
    count_parser.add_argument("-r", "--regions", required=True, type=str,
                              help="Genes or Peak File in BED, narrowPeak, or GTF Format")
    count_parser.add_argument("-b", "--barcodes", type=str, help="Two row TSV mapping barcode to cell-group/cluster")
    count_parser.add_argument("-o", "--output", type=str, help="Output Directory", default=str(Path.cwd()))

    count_parser.add_argument("--nofilt", action='store_true',
                              help="Skip step that filters BAM reads in regions of interest")  # TODO filtering options
    count_parser.add_argument("--keeptemps", type=str, nargs="?", const=0,
                              help="Keep intermediate files created during prefiltering. Outputs to directory if given, otherwise outputs ")

    analysis_parser = subparser.add_parser("analysis", parents=[parent_parser])
    analysis_parser.add_argument("counts", help="Count TSV output from count tool")
    analysis_parser.add_argument("--min", type=int, help="Minimum allele count for analysis", default=10)
    analysis_parser.add_argument("-m", "--model", type=str, choices=["single", "linear", "binomial"],
                                 help="Analysis Model", default="single")
    analysis_parser.add_argument("-o", "--output", type=str, help="Output Directory", default=str(Path.cwd()))

    args = parser.parse_args()

    print(args)

    return args


def run(args):
    if args.singlecell: # Single Cell Data
        print("Single Cell Analysis")
        if args.command == "count":
            print("Count Single-Cell")
            # parse_counting_sc(args.alignment, args.genotypes, args.regions, args.sample, args.barcodes, args.output, args.stype, nofilt=args.nofilt, temp_loc=args.keeptemps, features=args.features)
        
        elif args.command == "analysis":
            print("Analysis Single-Cell")
            # parse_analysis_sc(args.counts, args.min, args.model, args.output, args.stype, features=args.features)

    else: # Bulk processing
        print("Bulk Analysis")
        if args.command == "count":
            print("Count Bulk")
            # parse_counting(args.alignment, args.genotypes, args.regions, args.sample, args.output, args.stype, nofilt=args.nofilt, temp_loc=args.keeptemps, features=args.features)

        elif args.command == "analysis":
            print("Analysis Bulk")
            # parse_analysis(args.counts, args.min, args.model, args.output, args.stype, features=args.features)



def main():
    parent_parser = argparse.ArgumentParser(add_help=False)
    seq_type = parent_parser.add_mutually_exclusive_group() # Denote rna-seq vs atac-seq
    seq_type.add_argument("--rna", action='store_const', const="rna", dest="stype", help="Explicitly denote RNA-seq data (Otherwise infers from input)")
    seq_type.add_argument("--atac", action='store_const', const="atac", dest="stype", help="Explicitly denote ATAC-seq data (Otherwise infers from input)")

    parent_parser.add_argument("-sc", "--singlecell", action="store_true", help="Single Cell Option")

    parent_parser.add_argument("-ft", "--features", nargs="*", help="GTF features to analyze")

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

    count_parser.add_argument("--nofilt", action='store_true', help="Skip step that filters BAM reads in regions of interest") # TODO filtering options
    count_parser.add_argument("--keeptemps", type=str, nargs="?", const=0, help="Keep intermediate files created during prefiltering. Outputs to directory if given, otherwise outputs ")


    analysis_parser = subparser.add_parser("analysis", parents=[parent_parser])
    analysis_parser.add_argument("counts", help="Count TSV output from count tool")
    analysis_parser.add_argument("--min", type=int, help="Minimum allele count for analysis", default=10)
    analysis_parser.add_argument("-m", "--model", type=str, choices=["single", "linear", "binomial"], help="Analysis Model", default="single")
    analysis_parser.add_argument("-o", "--output", type=str, help="Output Directory", default=str(Path.cwd()))

    args = parser.parse_args()

    validate_out = validate_args(args)
    
    if isinstance(validate_out, str):
        parser.error(validate_out)
    else:
        args = validate_out
    
    print(args)
    
    if args.singlecell: # Single Cell Data
        print("Single Cell Analysis")
        if args.command == "count":
            parse_counting_sc(args.alignment, args.genotypes, args.regions, args.sample, args.barcodes, args.output, args.stype, nofilt=args.nofilt, temp_loc=args.keeptemps, features=args.features)
        
        elif args.command == "analysis":
            parse_analysis_sc(args.counts, args.min, args.model, args.output, args.stype, features=args.features)

    else: # Bulk processing
        print("Bulk Analysis")
        if args.command == "count":
            parse_counting(args.alignment, args.genotypes, args.regions, args.sample, args.output, args.stype, nofilt=args.nofilt, temp_loc=args.keeptemps, features=args.features)

        elif args.command == "analysis":
            parse_analysis(args.counts, args.min, args.model, args.output, args.stype, features=args.features)


if __name__ == '__main__':
    main()
