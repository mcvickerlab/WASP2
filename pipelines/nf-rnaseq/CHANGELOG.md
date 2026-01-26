# Changelog

All notable changes to the nf-rnaseq pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-01-25

### Added
- Initial release of WASP2 RNA-seq Allele-Specific Expression (ASE) pipeline
- WASP2 integration for mapping bias correction
- STAR aligner support with two-pass mapping
- Comprehensive samplesheet validation with edge case handling
- VCF index validation
- Allelic imbalance statistical testing with binomial test
- Skip analysis parameter (`--skip_analysis`) for optional imbalance testing
- Multiple output formats: TSV, Parquet, AnnData (H5AD), Zarr
- Lightweight real-data integration test suite
- WASP2 allelic analysis output validation tests
- nf-core compatible DSL2 module structure
- MultiQC integration for quality control reporting
- Support for Conda, Docker, and Singularity containers
