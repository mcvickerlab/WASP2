# Changelog

All notable changes to the nf-atacseq pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-01-25

### Added
- Initial release of WASP2 ATAC-seq Allelic Imbalance pipeline
- WASP2 integration for mapping bias correction in chromatin accessibility data
- BWA-MEM2 aligner support for ATAC-seq reads
- Peak calling with MACS2
- Full nf-core subworkflow pattern compliance
- Comprehensive meta.yml documentation for modules and subworkflows
- Validation test suite with edge case coverage
- Multiple output formats: TSV, BED, Parquet
- nf-core compatible DSL2 module structure
- MultiQC integration for quality control reporting
- Support for Conda, Docker, and Singularity containers
