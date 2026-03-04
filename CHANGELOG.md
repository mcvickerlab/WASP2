# Changelog

All notable changes to WASP2 will be documented in this file.

## [1.4.0] - 2026-02-25

### Fixed
- **14 verified bug fixes** across Python and Rust pipelines:
  - `not` → `~` on pandas Series for chromosome filtering (#50)
  - `samples[0]` truncation — pass full list instead of first element (#51)
  - `--phased`/`--region_col`/`--groupby` forwarded to Rust analysis (#52)
  - `AttributeError` on `is_gene_file` when no feature file provided (#53)
  - `NameError` on `json_dict` in mapping pipeline (#54)
  - Dead expression in `compare_ai.py` snp_cutoff calculation (#55)
  - VCF genotype Debug-format parsing fragility in Rust (#56)
  - `remap_names_path` checked in parallel fallback path (#57)
  - Removed `removed_moved` double-counting in mapping stats (#58)
  - BGZF reader fallback for standard gzip VCF files (#59)
  - uint8 → uint16 overflow in single-cell sparse matrices (#60)
  - Non-deterministic output order from parallel analysis (#62)
  - Phased check strict equality → issubset (#63)
  - `.all()` → `.any()` in GTF attribute detection (#64)
  - Phase detection reads 100 VCF records instead of 1 (#65)
- Linear model shape mismatch in analysis benchmarks
- `test_opt_prob` missing argument in test suite

### Added
- **pixi.toml** for one-command reproducible environments (`pixi install`)
- **Python 3.13** support (wheels built for 3.10–3.13)
- **Bioconda recipe** with maturin+Rust build support
- Mamba/conda installation documentation with miniforge link

### Changed
- Installation docs reordered: pixi (recommended) → pip → mamba/conda → Docker
- `release.yml` builds 16 wheels (4 Python × 4 platforms, up from 12)
- Version consistency checks cover Dockerfile and Singularity.def

## [1.3.0] - 2025-01-29

### Added
- **Nextflow pipeline ecosystem** with nf-core standards compliance
  - nf-atacseq pipeline for ATAC-seq allelic imbalance analysis
  - nf-rnaseq pipeline with validation and test suite
  - nf-scatac pipeline for single-cell ATAC-seq analysis
  - nf-outrider pipeline for WASP2+OUTRIDER integration
- **Distribution packaging** for PyPI, Bioconda, and Galaxy (#82)
- meta.yml documentation for Nextflow modules and subworkflows (#58)
- Validation test suites for all Nextflow pipelines (#51, #52, #54)

### Changed
- Nextflow modules now follow full nf-core subworkflow pattern compliance (#55, #60)
- Enhanced error handling in ATAC-seq Nextflow modules with warning logging and explicit error propagation
- Updated sample VCF test data files for better test coverage

### Fixed
- INDEL counting logic in Rust module (synced from WASP2-exp branch)
- Pandas and anndata version constraints for compatibility (#68)
- nf-core module robustness issues identified in PR review (parameter types, VCF index documentation)

## [1.2.0] - 2025-01-23

### Added
- **61× faster WASP filtering** via Rust optimization (validated r² > 0.99 vs GATK)
- INDEL support in variant processing
- bcftools and samtools added to environment.yml
- nf-test infrastructure for Nextflow modules
- Docker container support with ghcr.io publishing
- Security scanning workflow (pip-audit, cargo-audit, Bandit, Gitleaks, CodeQL)

### Fixed
- VCF→BED conversion now handles missing genotypes correctly
- CI maturin build fixed by using virtualenv
- Polars version constraint for stable API

### Changed
- Pinned pandas<2.0 and anndata<0.10 for compatibility
- Added ruff linting and pre-commit hooks for code quality
- Nextflow modules now use containerized WASP2

## [1.1.0] - 2024-11-24

### Added
- **Rust acceleration** for counting, mapping, and analysis modules (10-50x speedup)
- PyO3 bindings for seamless Python-Rust integration
- Multi-threaded BAM processing via `WASP2_RUST_THREADS` env var
- GitHub Pages documentation with PyData theme
- Validation scripts for parity testing

### Changed
- CLI now routes through Rust by default (no Python fallback for core operations)
- Updated to maturin-based build system for wheel packaging
- Modernized Sphinx docs with autodoc API generation

### Fixed
- Memory efficiency improvements in large BAM processing
- Consistent allele counting behavior across threads

## [1.0.0] - 2024-09-01

### Added
- Initial release
- Allele-specific read counting from BAM files
- WASP mapping bias correction algorithm
- Beta-binomial allelic imbalance analysis
- Single-cell allele counting support
- CLI tools: `wasp2-count`, `wasp2-map`, `wasp2-analyze`
