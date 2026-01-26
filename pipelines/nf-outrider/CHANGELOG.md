# Changelog

All notable changes to the nf-outrider pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-01-25

### Added
- Initial release of WASP2 + OUTRIDER pipeline for aberrant expression detection
- WASP2 integration for allele-specific expression analysis
- OUTRIDER algorithm for statistical outlier detection
- Mono-allelic expression (MAE) analysis
- Aberrant expression calling with configurable thresholds
- Multiple output formats: TSV, Parquet, AnnData (H5AD), Zarr
- meta.yml documentation for modules and subworkflows
- nf-core compatible DSL2 module structure
- MultiQC integration for quality control reporting
- Support for Conda, Docker, and Singularity containers
