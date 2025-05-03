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
    """
    A class to store and process file paths and parameters for WASP analysis.

    This class holds the input file paths (BAM and VCF), paired-end and phasing information,
    sample information, and output directories. It also generates default names for intermediary
    files required for downstream WASP analysis.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    is_paired : bool, optional
        Flag indicating whether the BAM file is paired-end. If None, it is determined automatically.
    samples : str or list of str, optional
        Sample information. If provided as a string, it is interpreted either as a file path
        (one sample per line) or as a comma-delimited list.
    is_phased : bool, optional
        Flag indicating whether the VCF file is phased. If None, it is determined by inspecting the VCF.
    out_dir : str, optional
        Output directory for final results. Defaults to the BAM file's directory if not provided.
    temp_loc : str, optional
        Directory for temporary files. Defaults to out_dir if not provided.

    Attributes
    ----------
    bam_file : str
        Path to the BAM file.
    vcf_file : str
        Path to the VCF file.
    is_paired : bool
        Whether the BAM file is paired-end.
    samples : list of str or None
        Processed list of sample names.
    is_phased : bool
        Whether the VCF file is phased.
    out_dir : str
        Output directory for results.
    temp_loc : str
        Directory for temporary files.
    vcf_prefix : str
        Prefix derived from the VCF filename.
    bam_prefix : str
        Prefix derived from the BAM filename.
    vcf_bed : str
        Path to the intermediary BED file converted from the VCF.
    remap_reads : str
        Path to a file listing read names to be remapped.
    intersect_file : str
        Path to the file containing intersections between VCF variants and regions.
    to_remap_bam : str
        Path to the BAM file containing reads to be remapped.
    keep_bam : str
        Path to the BAM file containing reads that are kept.
    remap_fq1 : str
        For paired-end data, path to the FASTQ file for read1.
    remap_fq2 : str or None
        For paired-end data, path to the FASTQ file for read2; None for single-end.
    """
    def __init__(
        self,
        bam_file: str,
        vcf_file: str,
        is_paired: Optional[bool] = None,
        samples: Optional[Union[str, List[str]]] = None,
        is_phased: Optional[bool] = None,
        out_dir: Optional[str] = None,
        temp_loc: Optional[str] = None
    ) -> None:
        # User input files
        self.bam_file: str = bam_file
        self.vcf_file: str = vcf_file
        self.is_paired: Optional[bool] = is_paired
        self.samples: Optional[Union[str, List[str]]] = samples
        self.is_phased: Optional[bool] = is_phased
        self.out_dir: Optional[str] = out_dir
        self.temp_loc: Optional[str] = temp_loc
        
        # Autoparse args
        if self.is_paired is None:
            with AlignmentFile(self.bam_file, "r") as bam:
                self.is_paired = next(bam.head(1)).is_paired
        
        # Process samples as list
        if self.samples is None:
            self.is_phased = False  # No phasing w/o sample
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
                samps_phased = [vcf_samps[s].phased for s in self.samples]  # assuming self.samples is a list
                if all(samps_phased):
                    self.is_phased = True
                else:
                    # TODO GOTTA WARN UNPHASED BAD
                    # TODO WARN SOME UNPHASED WHILE OTHERS PHASED
                    self.is_phased = False
        
        if self.out_dir is None:
            self.out_dir = str(Path(bam_file).parent)  # change to cwd?
        
        # TODO handle temp loc, maybe make default if temp not made?
        # Temporary workaround until figure out temp dir options
        if self.temp_loc is None:
            self.temp_loc = self.out_dir
        
        # Generate intermediate files
        # Maybe use easy default names if temp loc in use
        vcf_prefix = re.split(r'.vcf|.bcf', Path(self.vcf_file).name)[0]
        bam_prefix = Path(self.bam_file).name.rsplit(".bam")[0]
        
        self.vcf_prefix: str = vcf_prefix
        self.bam_prefix: str = bam_prefix
        
        self.vcf_bed: str = str(Path(self.temp_loc) / f"{vcf_prefix}.bed")
        self.remap_reads: str = str(Path(self.temp_loc) / f"{bam_prefix}_remap_reads.txt")
        self.intersect_file: str = str(Path(self.temp_loc) / f"{bam_prefix}_{vcf_prefix}_intersect.bed")
        
        self.to_remap_bam: str = str(Path(self.out_dir) / f"{bam_prefix}_to_remap.bam")
        self.keep_bam: str = str(Path(self.out_dir) / f"{bam_prefix}_keep.bam")
        
        # Relevant output reads
        if self.is_paired:
            self.remap_fq1: str = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles_r1.fq")
            self.remap_fq2: Optional[str] = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles_r2.fq")
        else:
            self.remap_fq1: str = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles.fq")
            self.remap_fq2 = None
    
    def write_data(self, out_file: Optional[str] = None) -> None:
        """
        Export relevant file data to a JSON file for downstream processing.

        This method exports the internal state of the WaspDataFiles object as a JSON file. The JSON
        file contains all relevant attributes needed for post-remapping analysis.

        Parameters
        ----------
        out_file : str, optional
            The output JSON file path. If not provided, defaults to "[BAM_PREFIX]_wasp_data_files.json"
            in the output directory.

        Returns
        -------
        None

        Side Effects
        ------------
        Writes a JSON file containing the object's __dict__.
        """
        if out_file is None:
            out_file = str(Path(self.out_dir) / f"{self.bam_prefix}_wasp_data_files.json")
        
        with open(out_file, "w") as json_out:
            json.dump(self.__dict__, json_out)
        
        print(f"File Data written to JSON...\n{out_file}")
