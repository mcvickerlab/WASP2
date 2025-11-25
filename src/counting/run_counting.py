import sys
import timeit
import re
import functools
import tempfile
import warnings

from pathlib import Path

# local imports
from .filter_variant_data import vcf_to_bed, intersect_vcf_region, parse_intersect_region, parse_intersect_region_new
from .parse_gene_data import make_gene_data, parse_intersect_genes, parse_intersect_genes_new
from .count_alleles import make_count_df

# Should I put this in separate file?
class WaspCountFiles:

    def __init__(self, bam_file, variant_file,
                 region_file=None, samples=None,
                 use_region_names=False,
                 out_file=None,
                 temp_loc=None,
                 precomputed_vcf_bed=None,
                 precomputed_intersect=None
                 ):

        # User input files
        self.bam_file = bam_file
        self.variant_file = variant_file
        self.region_file = region_file
        self.samples = samples
        self.use_region_names = use_region_names
        self.out_file = out_file
        self.temp_loc = temp_loc
        
        # gtf and gff specific
        self.is_gene_file = False # check if using gff3/gtf
        self.gtf_bed = None
        
        # Make sure samples turned into str list
        if isinstance(self.samples, str):
            
            # Check if sample file or comma delim string
            if Path(self.samples).is_file():
                
                with open(self.samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]
            
            else:
                self.samples = [s.strip() for s in self.samples.split(",")]
        
        
        # parse output?
        if self.out_file is None:
            self.out_file = str(Path.cwd() / "counts.tsv")
        
        
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

        # Filtered variant output (or precomputed)
        self.vcf_bed = precomputed_vcf_bed if precomputed_vcf_bed is not None else str(Path(self.temp_loc) / f"{variant_prefix}.bed")
        self.skip_vcf_to_bed = precomputed_vcf_bed is not None
        
        # Parse region file
        self.region_type = None # maybe use a boolean flag instead
        
        if self.region_file is not None:
            f_ext = "".join(Path(self.region_file).suffixes)
            
            if re.search(r'\.(.*Peak|bed)(?:\.gz)?$', f_ext, re.I):
                self.region_type = "regions"
                self.intersect_file = precomputed_intersect if precomputed_intersect is not None else str(Path(self.temp_loc) / f"{variant_prefix}_intersect_regions.bed")
                self.is_gene_file = False
            elif re.search(r'\.g[tf]f(?:\.gz)?$', f_ext, re.I):
                self.region_type = "genes"
                self.intersect_file = precomputed_intersect if precomputed_intersect is not None else str(Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.g[tf]f', Path(self.region_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_region_names = True # Use feature attributes as region names
            elif re.search(r'\.gff3(?:\.gz)?$', f_ext, re.I):
                self.region_type = "genes"
                self.intersect_file = precomputed_intersect if precomputed_intersect is not None else str(Path(self.temp_loc) / f"{variant_prefix}_intersect_genes.bed")
                self.is_gene_file = True
                gtf_prefix = re.split(r'.gff3', Path(self.region_file).name)[0]
                self.gtf_bed = str(Path(self.temp_loc) / f"{gtf_prefix}.bed")
                self.use_region_names = True # Use feature attributes as region names
            else:
                raise ValueError(f"Invalid region file type. Expected .bed, .gtf, or .gff3, got: {self.region_file}")

        else:
            # No region file: intersect file defaults to vcf_bed (or provided precomputed)
            self.intersect_file = precomputed_intersect if precomputed_intersect is not None else self.vcf_bed
        self.skip_intersect = precomputed_intersect is not None
        
        # TODO UPDATE THIS WHEN I ADD AUTOPARSERS
        if self.is_gene_file:
            
            # Possible edge case of vcf and gtf prefix conflict
            if self.vcf_bed == self.gtf_bed:
                self.gtf_bed = str(Path(self.temp_loc) / "genes.bed")
            

def tempdir_decorator(func):
    """Checks and makes tempdir for 
    run_count_variants()
    """
    
    @functools.wraps(func)
    def tempdir_wrapper(*args, **kwargs):
        
        if kwargs.get("temp_loc", None) is not None:
            return func(*args, **kwargs)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                kwargs["temp_loc"] = tmpdir
                return func(*args, **kwargs)
            
    return tempdir_wrapper


@tempdir_decorator
def run_count_variants(bam_file, variant_file,
                       region_file=None,
                       samples=None,
                       use_region_names=None,
                       out_file=None,
                       temp_loc=None,
                       gene_feature=None,
                       gene_attribute=None,
                       gene_parent=None,
                       use_rust=True,
                       precomputed_vcf_bed=None,
                       precomputed_intersect=None
                       ):


    # call the data class
    count_files = WaspCountFiles(bam_file, variant_file,
                                 region_file=region_file,
                                 samples=samples,
                                 use_region_names=use_region_names,
                                 out_file=out_file,
                                 temp_loc=temp_loc,
                                 precomputed_vcf_bed=precomputed_vcf_bed,
                                 precomputed_intersect=precomputed_intersect
                                )
    
    # print(*vars(count_files).items(), sep="\n") # For debugging
    with_gt = False
    if (count_files.samples is not None) and (len(count_files.samples) == 1):
        with_gt = True
        
        # temporarily disable for ASE
        # if not count_files.is_gene_file:
        #     with_gt = True
            
    
    # Create Intermediary Files
    if not count_files.skip_vcf_to_bed:
        vcf_to_bed(vcf_file=count_files.variant_file,
                   out_bed=count_files.vcf_bed,
                   samples=count_files.samples,
                   include_gt=with_gt
                  )
    
    
    # TODO PARSE GENE FEATURES AND ATTRIBUTES
    region_col_name = None # Defaults to 'region' as region name
    intersect_genes = False
    
    # region_files is valid to perform intersects
    if count_files.region_file is not None:
        
        # Check if we need to prepare genes for intersection
        if count_files.gtf_bed is not None:
            
            # TODO UPDATE THIS WHEN I ADD AUTOPARSERS AND VALIDATORS
            gene_data = make_gene_data(
                gene_file=count_files.region_file,
                out_bed=count_files.gtf_bed,
                feature=gene_feature,
                attribute=gene_attribute,
                parent_attribute=gene_parent
                )
            
            regions_to_intersect = count_files.gtf_bed
            region_col_name = gene_data.feature
            intersect_genes = True
        else:
            regions_to_intersect = count_files.region_file
            region_col_name = None # Defaults to 'region' as region name
        
        if not count_files.skip_intersect:
            intersect_vcf_region(vcf_file=count_files.vcf_bed,
                                 region_file=regions_to_intersect,
                                 out_file=count_files.intersect_file)
    

    # Create Variant Dataframe
    # TODO validate
    if intersect_genes:
        df = parse_intersect_genes_new(
            intersect_file=count_files.intersect_file,
            attribute=gene_data.attribute,
            parent_attribute=gene_data.parent_attribute)
    elif with_gt:
        df = parse_intersect_region_new(
            intersect_file=count_files.intersect_file,
            samples=["GT"],
            use_region_names=count_files.use_region_names,
            region_col=region_col_name
            )
    else:
        df = parse_intersect_region_new(
            intersect_file=count_files.intersect_file,
            samples=None,
            use_region_names=count_files.use_region_names,
            region_col=region_col_name
            )
        # df = parse_intersect_region(
        #     intersect_file=count_files.intersect_file,
        #     use_region_names=count_files.use_region_names,
        #     region_col=region_col_name)
    
    # Should I include a filt bam step???
    
    # Count
    count_df = make_count_df(bam_file=count_files.bam_file,
                             df=df,
                             use_rust=use_rust)
    
    # Write counts
    count_df.write_csv(count_files.out_file, include_header=True, separator="\t")
    
    # Should i return for use in analysis pipeline?
    # return count_df
