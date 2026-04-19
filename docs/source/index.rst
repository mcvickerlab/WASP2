WASP2: Allele-Specific Analysis
================================

.. image:: https://img.shields.io/pypi/v/wasp2
   :target: https://pypi.org/project/wasp2/
   :alt: PyPI

.. image:: https://github.com/mcvickerlab/WASP2/workflows/WASP2%20Tests/badge.svg
   :target: https://github.com/mcvickerlab/WASP2/actions
   :alt: Tests

WASP2 is a comprehensive suite of tools for unbiased allele-specific analysis of next-generation sequencing data. It addresses reference bias in read mapping and provides statistical methods for detecting allelic imbalance.

Features
--------

* **Unbiased Mapping**: WASP algorithm for correcting reference bias
  (van de Geijn et al. 2015, Nat Methods 10.1038/nmeth.3582)
* **Allele Counting**: Count allele-specific reads from BAM files
* **Statistical Analysis**: Beta-binomial likelihood-ratio tests for allelic
  imbalance
* **Single-Cell Support**: Specialized tools for single-cell RNA-seq and
  scATAC-seq
* **Rust-Accelerated**: Rust-backed BAM filtering and counting for large
  cohorts

Quick Start
-----------

Install via pip:

.. code-block:: bash

   pip install wasp2

Count alleles from a BAM file:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf

Analyze allelic imbalance:

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv

Documentation
-------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart
   choosing_workflow
   faq

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/counting
   user_guide/mapping
   user_guide/analysis
   user_guide/single_cell

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/quickstart_counting
   tutorials/bulk_workflow
   tutorials/scrna_seq
   tutorials/scatac_workflow
   tutorials/comparative_imbalance

.. toctree::
   :maxdepth: 2
   :caption: Statistical Methods

   methods/index
   methods/counting_algorithm
   methods/mapping_filter
   methods/statistical_models
   methods/dispersion_estimation
   methods/fdr_correction

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/counting
   api/mapping
   api/analysis

.. toctree::
   :maxdepth: 1
   :caption: Development

   development
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
