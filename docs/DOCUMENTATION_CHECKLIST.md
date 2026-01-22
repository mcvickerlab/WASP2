# WASP2 Documentation Implementation Checklist

Track progress on documentation improvements. Mark items as complete with [x].

## Phase 1: Quick Wins (1-2 weeks)

### README Enhancements
- [ ] Add enhanced badge section (CI, coverage, downloads, conda)
- [ ] Move Quick Start section before Installation
- [ ] Add Feature Highlights section with clear hierarchy
- [ ] Create Installation Options matrix (PyPI, conda, source, codespaces)
- [ ] Add Citation section with BibTeX
- [ ] Add Comparison Table (vs GATK, phASER, MBASED)
- [ ] Add Learning Path section linking to tutorials
- [ ] Test all README code blocks for accuracy

### Quick Reference Materials
- [ ] Create CHEATSHEET.md with common commands
- [ ] Add one-liner examples directory (examples/README.md)
- [ ] Create example shell scripts (basic_rnaseq.sh, basic_atacseq.sh)
- [ ] Add small test dataset for tutorials

### FAQ Section
- [ ] Create docs/source/faq.rst
- [ ] Add 10-15 most common questions
- [ ] Include troubleshooting Q&A
- [ ] Link from main documentation index

### Shell Completion
- [ ] Generate bash completion script
- [ ] Generate zsh completion script
- [ ] Generate fish completion script
- [ ] Add installation instructions to README
- [ ] Test completion scripts on each shell

---

## Phase 2: Core Documentation (2-3 weeks)

### Tutorial Series

#### Tutorial 0: Concepts
- [ ] Create docs/tutorials/00_concepts.md
- [ ] Explain allelic imbalance with examples
- [ ] Describe reference bias problem
- [ ] Illustrate WASP solution with diagram
- [ ] Add decision tree for when to use each module

#### Tutorial 1: Quick Start (5 min)
- [ ] Create docs/tutorials/01_quickstart.md
- [ ] Prepare small test dataset (~50 MB)
- [ ] Write 5-minute end-to-end example
- [ ] Test timing on fresh system
- [ ] Add expected outputs

#### Tutorial 2: Installation Guide
- [ ] Create docs/tutorials/02_installation_guide.md
- [ ] Cover all installation methods
- [ ] Add platform-specific instructions (Linux, macOS, Windows/WSL)
- [ ] Include troubleshooting common install issues
- [ ] Verify each installation method

#### Tutorial 3: Basic Workflow (30 min)
- [ ] Create docs/tutorials/03_basic_workflow.md
- [ ] Cover complete pipeline (QC → WASP → Count → Analyze)
- [ ] Add pipeline diagram
- [ ] Include interpretation section
- [ ] Add quality control checks

#### Tutorial 4: RNA-seq ASE (45 min)
- [ ] Create docs/tutorials/04_rnaseq_ase.md
- [ ] Use realistic dataset (GM12878 or similar)
- [ ] Cover gene-level analysis
- [ ] Include visualization examples
- [ ] Add validation against known imprinted genes

#### Tutorial 5: ATAC-seq ASE (45 min)
- [ ] Create docs/tutorials/05_atac_ase.md
- [ ] Cover peak calling integration
- [ ] Explain differences from RNA-seq
- [ ] Include TF motif enrichment section
- [ ] Add caQTL interpretation

#### Tutorial 6: Single-Cell (60 min)
- [ ] Create docs/tutorials/06_single_cell.md
- [ ] Cover 10x Genomics workflow
- [ ] Explain cell-type-specific analysis
- [ ] Include differential AI section
- [ ] Add visualization in Python (scanpy)

#### Tutorial 7: Advanced Options
- [ ] Create docs/tutorials/07_advanced_options.md
- [ ] Cover all command-line options
- [ ] Explain parameter tuning
- [ ] Include use case examples

#### Tutorial 8: Troubleshooting (reference)
- [ ] Create docs/tutorials/08_troubleshooting.md
- [ ] Organize by module (count, map, analyze)
- [ ] Add diagnostic commands for each issue
- [ ] Include error message reference table
- [ ] Add decision trees for common problems

#### Tutorial 9: Performance Tuning
- [ ] Create docs/tutorials/09_performance_tuning.md
- [ ] Benchmark different variant formats
- [ ] Explain threading and parallelization
- [ ] Cover memory optimization strategies
- [ ] Add HPC/cloud computing examples

### Enhanced CLI Help

#### Count Module
- [ ] Enhance count-variants help text with examples
- [ ] Enhance count-variants-sc help text
- [ ] Add output format descriptions
- [ ] Include performance tips in help

#### Map Module
- [ ] Enhance make-reads help text
- [ ] Enhance filter-remapped help text
- [ ] Add workflow diagram reference
- [ ] Include parameter recommendations

#### Analysis Module
- [ ] Enhance find-imbalance help text
- [ ] Enhance find-imbalance-sc help text
- [ ] Enhance compare-imbalance help text
- [ ] Add interpretation guidance

### CLI Reference Documentation
- [ ] Create docs/source/cli/index.rst
- [ ] Create docs/source/cli/wasp2_count.rst (complete reference)
- [ ] Create docs/source/cli/wasp2_map.rst
- [ ] Create docs/source/cli/wasp2_analyze.rst
- [ ] Add examples section to each
- [ ] Link from main documentation index

---

## Phase 3: Advanced Documentation (2-3 weeks)

### Man Pages

#### Main Man Pages
- [ ] Create man/man1/wasp2.1 (overview)
- [ ] Create man/man1/wasp2-count.1
- [ ] Create man/man1/wasp2-map.1
- [ ] Create man/man1/wasp2-analyze.1

#### Subcommand Man Pages
- [ ] Create man/man1/wasp2-count-variants.1
- [ ] Create man/man1/wasp2-count-variants-sc.1
- [ ] Create man/man1/wasp2-map-make-reads.1
- [ ] Create man/man1/wasp2-map-filter-remapped.1
- [ ] Create man/man1/wasp2-analyze-find-imbalance.1
- [ ] Create man/man1/wasp2-analyze-find-imbalance-sc.1
- [ ] Create man/man1/wasp2-analyze-compare-imbalance.1

#### Man Page Installation
- [ ] Add man pages to pyproject.toml data_files
- [ ] Test man page installation
- [ ] Verify man page formatting (groff)
- [ ] Test on different systems

### API Documentation (Comprehensive Docstrings)

#### Counting Module
- [ ] Add/enhance module docstring (counting/__init__.py)
- [ ] Enhance run_count_variants docstring
- [ ] Enhance run_count_variants_sc docstring
- [ ] Enhance WaspCountFiles docstring
- [ ] Add docstrings to all helper functions
- [ ] Run doctest on all examples

#### Mapping Module
- [ ] Add/enhance module docstring (mapping/__init__.py)
- [ ] Enhance run_make_remap_reads docstring
- [ ] Enhance run_wasp_filt docstring
- [ ] Add docstrings to all helper functions
- [ ] Run doctest on all examples

#### Analysis Module
- [ ] Add/enhance module docstring (analysis/__init__.py)
- [ ] Enhance run_ai_analysis docstring
- [ ] Enhance run_ai_analysis_sc docstring
- [ ] Enhance run_ai_comparison docstring
- [ ] Add docstrings to all statistical functions
- [ ] Run doctest on all examples

#### I/O Module
- [ ] Create comprehensive docstrings for VariantSource
- [ ] Document VCFSource, CyVCF2Source, PGENSource
- [ ] Add examples for each variant format
- [ ] Document performance characteristics

### Jupyter Notebook Examples
- [ ] Create examples/notebooks/basic_analysis.ipynb
- [ ] Create examples/notebooks/rnaseq_workflow.ipynb
- [ ] Create examples/notebooks/atacseq_workflow.ipynb
- [ ] Create examples/notebooks/visualization.ipynb
- [ ] Create examples/notebooks/single_cell_analysis.ipynb
- [ ] Test all notebooks execute without errors
- [ ] Add to documentation with nbsphinx

### Integration Guides
- [ ] Create how_to/integrate_with_nextflow.md
- [ ] Create how_to/integrate_with_snakemake.md
- [ ] Create how_to/integrate_with_cwl.md
- [ ] Create how_to/batch_processing.md
- [ ] Create how_to/cloud_deployment.md

---

## Phase 4: Polish (1 week)

### Visual Elements
- [ ] Create WASP algorithm diagram (SVG or PNG)
- [ ] Create pipeline flowchart
- [ ] Create decision tree for module selection
- [ ] Add before/after mapping bias illustration
- [ ] Create output format visual examples

### Enhanced Sphinx Documentation

#### Structure
- [ ] Create how_to/ directory and index
- [ ] Create explanations/ directory and index
- [ ] Create data_formats/ directory and index
- [ ] Reorganize existing pages to fit Divio structure
- [ ] Update navigation and cross-links

#### New Pages
- [ ] Create explanations/allelic_imbalance.rst
- [ ] Create explanations/reference_bias.rst
- [ ] Create explanations/wasp_algorithm.rst
- [ ] Create explanations/statistical_models.rst
- [ ] Create data_formats/input_formats.rst
- [ ] Create data_formats/output_formats.rst
- [ ] Create data_formats/variant_formats.rst
- [ ] Create how_to/interpret_results.rst

#### Enhancements
- [ ] Add sphinx-design cards to index page
- [ ] Add sphinx-tabs for format comparisons
- [ ] Add sphinx-copybutton configuration
- [ ] Enable myst_parser for Markdown support
- [ ] Add version switcher (if using RTD)

### Documentation Testing
- [ ] Set up documentation build in CI
- [ ] Add linkcheck to CI pipeline
- [ ] Add spell checking (optional)
- [ ] Test documentation builds on different Python versions
- [ ] Verify all code examples execute
- [ ] Run doctest on all docstrings

### Video Tutorials (Optional)
- [ ] Record 5-minute quick start screencast
- [ ] Record RNA-seq workflow walkthrough
- [ ] Record single-cell analysis demo
- [ ] Upload to YouTube
- [ ] Embed in documentation

---

## Ongoing Maintenance

### Version Management
- [ ] Set up Read the Docs with version switching
- [ ] Configure .readthedocs.yml
- [ ] Tag documentation versions with releases
- [ ] Maintain CHANGELOG.md
- [ ] Update docs/source/changelog.rst

### Quality Metrics
- [ ] Track docstring coverage (pydocstyle or interrogate)
- [ ] Monitor broken links (weekly check)
- [ ] Review GitHub issues tagged "documentation"
- [ ] Track most-searched terms (if analytics enabled)
- [ ] Collect user feedback

### Updates
- [ ] Update documentation with each release
- [ ] Keep performance benchmarks current
- [ ] Add new examples as features are added
- [ ] Refresh screenshots and outputs
- [ ] Review and update FAQ based on issues

---

## Priority Matrix

### High Priority (Do First)
1. Enhanced README (immediate value)
2. Quick Start tutorial (user onboarding)
3. FAQ section (reduce support burden)
4. Enhanced CLI help (daily use)
5. Basic workflow tutorial (complete pipeline)

### Medium Priority (Do Second)
1. Man pages (professional polish)
2. Comprehensive docstrings (API users)
3. RNA-seq and ATAC-seq tutorials (common workflows)
4. Troubleshooting guide (reduce support time)
5. Performance tuning guide (power users)

### Lower Priority (Nice to Have)
1. Video tutorials (multimedia learners)
2. Jupyter notebooks (interactive examples)
3. Pipeline integration guides (advanced users)
4. Additional visual diagrams (visual learners)
5. Translation (if international audience)

---

## Success Metrics

Track these to measure documentation effectiveness:

- [ ] Reduced "documentation" tagged issues
- [ ] Increased PyPI downloads after improvements
- [ ] Positive user feedback on tutorials
- [ ] Decreased response time on support questions
- [ ] Higher stars/forks on GitHub
- [ ] Citations in papers

---

## Resources Needed

### Tools
- [ ] Sphinx and extensions installed
- [ ] Documentation build environment
- [ ] Screen recording software (for videos)
- [ ] Diagram creation tool (draw.io, Inkscape, etc.)

### Data
- [ ] Test datasets for tutorials (<100 MB each)
- [ ] Example outputs for all commands
- [ ] Benchmark results for performance docs

### Time Estimates
- Phase 1 (Quick Wins): 10-15 hours
- Phase 2 (Core Docs): 30-40 hours
- Phase 3 (Advanced): 25-35 hours
- Phase 4 (Polish): 10-15 hours
- **Total**: 75-105 hours (2-3 months part-time)

---

## Notes

- Start with Phase 1 for immediate value
- Prioritize based on user feedback and common questions
- Iterate on tutorials with user testing
- Keep documentation version-controlled with code
- Update docs with every significant code change

**Last Updated**: 2025-01-22
