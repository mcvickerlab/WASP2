from pathlib import Path
import tempfile
import re
import json
from typing import Optional, Union, List

import pysam
from pysam import VariantFile
from pysam.libcalignmentfile import AlignmentFile


# TODO, GOTTA INCLUDE ALL POSSIBLE DATA COMBOS
class WaspDataFiles:
    """Manage file paths and auto-detection for WASP mapping pipeline."""

    def __init__(
        self,
        bam_file: Union[str, Path],
        vcf_file: Union[str, Path],
        is_paired: Optional[bool] = None,
        samples: Optional[Union[str, List[str]]] = None,
        is_phased: Optional[bool] = None,
        out_dir: Optional[Union[str, Path]] = None,
        temp_loc: Optional[Union[str, Path]] = None
    ) -> None:
        
        # User input files
        self.bam_file = bam_file
        self.vcf_file = vcf_file
        self.is_paired = is_paired
        self.samples = samples
        self.is_phased = is_phased
        self.out_dir = out_dir
        self.temp_loc = temp_loc
        
        
        # Autoparse args
        if self.is_paired is None:
            with AlignmentFile(self.bam_file, "r") as bam:
                self.is_paired = next(bam.head(1)).is_paired
        
        
        # Process samples as list
        if self.samples is None:
            self.is_phased = False # No phasing w/o sample
        elif isinstance(self.samples, str):
            
            # Check if sample file or comma delim string
            if Path(self.samples).is_file():
                
                with open(self.samples) as sample_file:
                    self.samples = [l.strip() for l in sample_file]
            
            else:
                self.samples = [s.strip() for s in self.samples.split(",")]
                # self.samples = self.samples.split(",") # should i strip spaces?
        
        # Check if VCF is phased
        if self.is_phased is None:
            # TODO GOTTA FIX THIS TO CHECK IF PHASED
            
            with VariantFile(self.vcf_file, "r") as vcf:
                vcf_samps = next(vcf.fetch()).samples
                samps_phased = [vcf_samps[s].phased for s in self.samples]
                
                if all(samps_phased):
                    self.is_phased = True
                else:
                    # TODO GOTTA WARN UNPHASED BAD
                    # TODO WARN SOME UNPHASED WHILE OTHERS PHASED
                    self.is_phased = False
        
        if self.out_dir is None:
            self.out_dir = Path(bam_file).parent # change to cwd?
        
        # TODO handle temp loc, maybe make default if temp not made?
        # Temporary workaround until figure out temp dir options
        if self.temp_loc is None:
            self.temp_loc = self.out_dir
        
        # Generate intermediate files
        # Maybe use easy defalt names if temp loc in use
        
        vcf_prefix = re.split(r'.vcf|.bcf', Path(self.vcf_file).name)[0]
        bam_prefix = Path(self.bam_file).name.rsplit(".bam")[0]
        
        self.vcf_prefix = vcf_prefix
        self.bam_prefix = bam_prefix
        
        self.vcf_bed = str(Path(self.temp_loc) / f"{vcf_prefix}.bed")
        self.remap_reads = str(Path(self.temp_loc) / f"{bam_prefix}_remap_reads.txt")
        self.intersect_file = str(Path(self.temp_loc) / f"{bam_prefix}_{vcf_prefix}_intersect.bed")
        
        self.to_remap_bam = str(Path(self.out_dir) / f"{bam_prefix}_to_remap.bam")
        self.keep_bam = str(Path(self.out_dir) / f"{bam_prefix}_keep.bam")
        
        # Relevant output reads
        if self.is_paired:
            self.remap_fq1 = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles_r1.fq")
            self.remap_fq2 = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles_r2.fq")
        else:
            self.remap_fq1 = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles.fq")
            self.remap_fq2 = None
    
    def write_data(self, out_file: Optional[Union[str, Path]] = None) -> None:
        """Export Relevant Files to JSON
        Used for parsing post remapping step easily

        :param out_file: name for output file if not using default
        :type out_file: str, optional
        """
        
        if out_file is None:
            out_file = str(Path(self.out_dir) / f"{self.bam_prefix}_wasp_data_files.json")
        
        with open(out_file, "w") as json_out:
            json.dump(self.__dict__, json_out)
        
        print(f"File Data written to JSON...\n{out_file}")
