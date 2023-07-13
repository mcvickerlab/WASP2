import sys
import timeit
import re
import functools
import tempfile
import warnings

from pathlib import Path

# local imports
from filter_variant_data import vcf_to_bed, intersect_vcf_region, parse_intersect_region
from count_alleles import make_count_df

# Should I put this in separate file?
class WaspCountFiles:

    def __init__(self, bam_file, vcf_file,
                 region_file=None, samples=None,
                 use_region_names=False,
                 out_file=None,
                 out_dir=None,
                 temp_loc=None):
        
        # User input files
        self.bam_file = bam_file
        self.vcf_file = vcf_file
        self.region_file = region_file
        self.samples = samples
        self.use_region_names = use_region_names
        self.out_file = out_file
        self.out_dir = out_dir  # should i replace this with out file???
        self.temp_loc = temp_loc
        
        
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
        
        # Parse vcf and intersect
        vcf_prefix = re.split(r'.vcf|.bcf', Path(self.vcf_file).name)[0]
        self.vcf_prefix = vcf_prefix
        
        # Filtered vcf output
        self.vcf_bed = str(Path(self.temp_loc) / f"{vcf_prefix}.bed")
        
        
        # Parse region file
        self.region_type = None # maybe use a boolean flag instead
        # self.intersect_file = None
        
        if self.region_file is not None:
            f_ext = "".join(Path(self.region_file).suffixes)
            
            if re.search(r'\.(.*Peak|bed)(\.gz)?$', f_ext, re.I):
                self.region_type = "regions"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_regions.bed")
                # file_stype = "atac"
            elif re.search(r'\.(gtf)(\.gz)?$', f_ext, re.I):
                # TODO gotta implement this
                print("NOT IMPLEMENTED YET")
                self.region_type = "genes"
                self.intersect_file = str(Path(self.temp_loc) / f"{vcf_prefix}_intersect_genes.bed")
                # file_stype = "rna" # TODO MAKE THIS WORK/gff too?
            else:
                self.region_file = None
                print("invalid ftype") # Make this raise an error later

        else:
            self.intersect_file = self.vcf_bed


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
def run_count_variants(bam_file, vcf_file,
                       region_file=None,
                       samples=None,
                       use_region_names=None,
                       out_file=None,
                       temp_loc=None):
    
    
    # call the data class
    count_files = WaspCountFiles(bam_file, vcf_file,
                                 region_file=region_file,
                                 samples=samples,
                                 use_region_names=use_region_names,
                                 out_file=out_file,
                                 temp_loc=temp_loc
                                )
    
    # Create Intermediary Files
    vcf_to_bed(vcf_file=count_files.vcf_file,
               out_bed=count_files.vcf_bed,
               samples=count_files.samples,
               drop_gt=True
              )
    
    
    # TODO GOTTA HANDLE GENES
    
    # region_files is valid to perform intersects
    if count_files.region_file is not None:
    
        intersect_vcf_region(vcf_file=count_files.vcf_bed,
                             region_file=count_files.region_file,
                             out_file=count_files.intersect_file)
    

    # Create Variant Dataframe
    df = parse_intersect_region(
        intersect_file=count_files.intersect_file,
        use_region_names = count_files.use_region_names
    )
    
    
    # Should I include a filt bam step???
    
    # Count
    count_df = make_count_df(bam_file=count_files.bam_file,
                             df=df)
    
    # Write counts
    count_df.write_csv(count_files.out_file, has_header=True, separator="\t")
    
    # Should i return for use in analysis pipeline?
    # return count_df