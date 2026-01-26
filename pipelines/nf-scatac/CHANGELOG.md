# Changelog

All notable changes to the nf-scatac pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-01-25

### Added
- Initial release of WASP2 Single-Cell ATAC-seq Allelic Imbalance pipeline
- True allele-specific counting with ref/alt/hap1/hap2 layers (BAM input mode)
- Overlap counting support for fragments input mode
- 10x Genomics CellRanger ATAC output support
- Cell barcode filtering with customizable thresholds
- Peak region filtering
- Per-cell allele counts with AnnData (H5AD) output
- AnnData layers: X (total), ref, alt, hap1, hap2 (when BAM provided)
- Zarr output format for GenVarLoader integration
- Pseudo-bulk aggregation for statistical power
- Cell QC metrics reporting
- Allelic imbalance statistical analysis
- nf-core subworkflow pattern compliance
- Validation and test suite
- Integration support for ArchR/Signac via scverse ecosystem
- Support for Conda, Docker, and Singularity containers
