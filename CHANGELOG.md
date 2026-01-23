# Changelog

All notable changes to WASP2 will be documented in this file.

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
