"""File path management for WASP mapping pipeline.

Provides the WaspDataFiles class for managing input/output paths
and auto-detecting file properties.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

from pysam import VariantFile
from pysam.libcalignmentfile import AlignmentFile

logger = logging.getLogger(__name__)


class WaspDataFiles:
    """Manage file paths and auto-detection for WASP mapping pipeline."""

    def __init__(
        self,
        bam_file: str | Path,
        variant_file: str | Path,
        is_paired: bool | None = None,
        samples: str | list[str] | None = None,
        is_phased: bool | None = None,
        out_dir: str | Path | None = None,
        temp_loc: str | Path | None = None,
    ) -> None:
        # User input files
        self.bam_file = bam_file
        self.variant_file = variant_file
        self.is_paired = is_paired
        self.samples = samples
        self.is_phased = is_phased
        self.out_dir = out_dir
        self.temp_loc = temp_loc

        # Autoparse args
        if self.is_paired is None:
            with AlignmentFile(str(self.bam_file), "r") as bam:
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

        # At this point, self.samples is normalized to Optional[List[str]]

        # Check if variant file is phased (only works for VCF/BCF, not PGEN)
        if self.is_phased is None and self.samples is not None:
            # TODO GOTTA FIX THIS TO CHECK IF PHASED
            # Note: This only works for VCF/BCF files, PGEN doesn't store phase in the same way
            variant_path = Path(self.variant_file)
            suffix = variant_path.suffix.lower()
            if suffix in (".vcf", ".bcf") or str(variant_path).lower().endswith(".vcf.gz"):
                with VariantFile(str(self.variant_file), "r") as vcf:
                    n_phased = 0
                    n_checked = 0
                    for rec in vcf.fetch():
                        n_checked += 1
                        if all(rec.samples[s].phased for s in self.samples):
                            n_phased += 1
                        if n_checked >= 100:
                            break

                    if n_checked > 0 and n_phased > n_checked // 2:
                        self.is_phased = True
                    else:
                        # TODO GOTTA WARN UNPHASED BAD
                        # TODO WARN SOME UNPHASED WHILE OTHERS PHASED
                        self.is_phased = False
            else:
                # PGEN format - assume phased (user should specify if not)
                self.is_phased = True

        if self.out_dir is None:
            self.out_dir = Path(bam_file).parent  # change to cwd?

        # TODO handle temp loc, maybe make default if temp not made?
        # Temporary workaround until figure out temp dir options
        if self.temp_loc is None:
            self.temp_loc = self.out_dir

        # Generate intermediate files
        # Maybe use easy defalt names if temp loc in use

        # Handle different variant file extensions for prefix extraction
        variant_name = Path(self.variant_file).name
        if variant_name.endswith(".vcf.gz"):
            variant_prefix = variant_name[:-7]  # Remove .vcf.gz
        elif variant_name.endswith(".pgen"):
            variant_prefix = variant_name[:-5]  # Remove .pgen
        else:
            variant_prefix = re.split(r"\.vcf|\.bcf", variant_name)[0]
        bam_prefix = Path(self.bam_file).name.rsplit(".bam")[0]

        self.variant_prefix = variant_prefix
        self.bam_prefix = bam_prefix

        self.vcf_bed = str(Path(self.temp_loc) / f"{variant_prefix}.bed")
        self.remap_reads = str(Path(self.temp_loc) / f"{bam_prefix}_remap_reads.txt")
        self.intersect_file = str(
            Path(self.temp_loc) / f"{bam_prefix}_{variant_prefix}_intersect.bed"
        )

        self.to_remap_bam = str(Path(self.out_dir) / f"{bam_prefix}_to_remap.bam")
        self.keep_bam = str(Path(self.out_dir) / f"{bam_prefix}_keep.bam")

        # Relevant output reads
        if self.is_paired:
            self.remap_fq1 = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles_r1.fq")
            self.remap_fq2: str | None = str(
                Path(self.out_dir) / f"{bam_prefix}_swapped_alleles_r2.fq"
            )
        else:
            self.remap_fq1 = str(Path(self.out_dir) / f"{bam_prefix}_swapped_alleles.fq")
            self.remap_fq2 = None

    def write_data(self, out_file: str | Path | None = None) -> None:
        """Export Relevant Files to JSON
        Used for parsing post remapping step easily

        :param out_file: name for output file if not using default
        :type out_file: str, optional
        """
        if out_file is None:
            out_file = str(Path(str(self.out_dir)) / f"{self.bam_prefix}_wasp_data_files.json")

        with open(out_file, "w") as json_out:
            json.dump(self.__dict__, json_out)

        logger.info("File data written to JSON: %s", out_file)
