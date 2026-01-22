import sys
import timeit
import re
import functools
import tempfile
import warnings

from pathlib import Path

import anndata as ad


# local imports
from .filter_variant_data import vcf_to_bed, intersect_vcf_region, parse_intersect_region_new
from .run_counting import tempdir_decorator
from .count_alleles_sc import make_count_matrix


class WaspCountSC:

    def __init__(self, bam_file,
                 variant_file,
                 barcode_file,
                 feature_file,
                 samples=None,
                 use_region_names=False,
                 out_file=None,
                 temp_loc=None
                 ):

        # TODO: ALSO ACCEPT .h5

        # User input files
        self.bam_file = bam_file
        self.variant_file = variant_file
        self.barcode_file = barcode_file # Maybe could be optional?
        
        self.feature_file = feature_file
        self.samples = samples
        self.use_region_names = use_region_names
        self.out_file = out_file
        self.temp_loc = temp_loc
        
        # Optional inputs and outputs?
        # output_sparse_mtx = None
        # SNP OUTPUT?!?!?

        
        # Make sure samples turned into str list
        # Ideally single sample for single cell
        if isinstance(self.samples, str):
            
            # Check if sample file or comma delim string
            if Path(self.samples).is_file():
                
                with open(self.samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]
            
            else:
                self.samples = [s.strip() for s in self.samples.split(",")]
        
        # parse output?
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "allele_counts.h5ad")
        
        
        # Failsafe if decorator doesnt create temp_loc
        if self.temp_loc is None:
            self.temp_loc = str(Path.cwd())
        
        
        # Parse variant file prefix (handle VCF, BCF, PGEN)
        variant_name = Path(self.variant_file).name
        if variant_name.endswith('.vcf.gz'):
            variant_prefix = variant_name[:-7]  # Remove .vcf.gz
        elif variant_name.endswith('.pgen'):
            variant_prefix = variant_name[:-5]  # Remove .pgen
        else:
            variant_prefix = re.split(r'\.vcf|\.bcf', variant_name)[0]
        self.variant_prefix = variant_prefix

        # Filtered variant output
        self.vcf_bed = str(Path(self.temp_loc) / f"{variant_prefix}.bed")
        
        # Parse feature file
        self.feature_type = None # maybe use a boolean flag instead
        
        if self.feature_file is not None:
            
            f_ext = "".join(Path(self.feature_file).suffixes)
            
            if re.search(r'\.(.*Peak|bed)(?:\.gz)?$', f_ext, re.I):
                self.feature_type = "regions"
                self.intersect_file = str(Path(self.temp_loc) / f"{variant_prefix}_intersect_regions.bed")
                self.is_gene_file = False
            elif re.search(r'\.g[tf]f(?:\.gz)?$', f_ext, re.I):
                self.feature_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.g[tf]f', Path(self.feature_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_feature_names = True # Use feature attributes as region names
            elif re.search(r'\.gff3(?:\.gz)?$', f_ext, re.I):
                self.feature_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.gff3', Path(self.feature_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_feature_names = True # Use feature attributes as feature names
            else:
                raise ValueError(f"Invalid feature file type. Expected .bed, .gtf, or .gff3, got: {self.feature_file}")

        else:
            self.intersect_file = self.vcf_bed
        
        # TODO UPDATE THIS WHEN I ADD AUTOPARSERS
        if self.is_gene_file:
            
            # Possible edge case of vcf and gtf prefix conflict
            if self.vcf_bed == self.gtf_bed:
                self.gtf_bed = str(Path(self.temp_loc) / "genes.bed")



@tempdir_decorator
def run_count_variants_sc(bam_file, variant_file,
                          barcode_file,
                          feature_file=None,
                          samples=None,
                          use_region_names=False,
                          out_file=None,
                          temp_loc=None,
                         ):

    # Stores file names
    count_files = WaspCountSC(bam_file, variant_file,
                              barcode_file=barcode_file,
                              feature_file=feature_file,
                              samples=samples,
                              use_region_names=use_region_names,
                              out_file=out_file,
                              temp_loc=temp_loc
                             )

    # Create intermediary files
    # Maybe change include_gt based on preparse?
    vcf_to_bed(vcf_file=count_files.variant_file,
               out_bed=count_files.vcf_bed,
               samples=count_files.samples,
               include_gt=True
              )

    intersect_vcf_region(vcf_file=count_files.vcf_bed,
                         region_file=count_files.feature_file,
                         out_file=count_files.intersect_file)


    # TODO: handle use_region_names better
    df = parse_intersect_region_new(
        intersect_file=count_files.intersect_file,
        samples=count_files.samples,
        use_region_names=use_region_names,
        region_col=None
    )
    

    # TODO: handle case where barcode file contains multiple columns
    with open(count_files.barcode_file, "r") as file:
        bc_dict = {
            line.rstrip():i for i, line in enumerate(file)}
    
    # Generate Output
    adata = make_count_matrix(bam_file=count_files.bam_file,
                              df=df, bc_dict=bc_dict,
                              include_samples=count_files.samples
                             )
    
    # Write outputs
    adata.write_h5ad(count_files.out_file)
    # TODO: include output options, (ie MTX, dense?)

    