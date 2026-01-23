Changelog
=========

Version 1.0.0 (2025-11-17)
--------------------------

Initial Release
~~~~~~~~~~~~~~~

**Features:**

* Complete type hint coverage (24 files, 5,500 lines)
* PyPI package available (pip install wasp2)
* CI/CD pipeline with GitHub Actions
* Pre-commit hooks for code quality
* Comprehensive documentation on ReadTheDocs

**Modules:**

* **Counting**: Allele-specific read counting from BAM files
* **Mapping**: WASP algorithm for unbiased read remapping
* **Analysis**: Statistical detection of allelic imbalance

**Type Hints:**

* TH-1: Counting module (7 files)
* TH-2: Analysis module (10 files)
* TH-3: Mapping module (7 files)

**Testing:**

* Regression tests (memory, performance)
* Full pipeline validation with real genomic data
* All tests passing in CI

**Documentation:**

* API documentation auto-generated from type hints
* User guides for each module
* Installation and quickstart guides
* Development and contributing guides
