WASP2: Allele-Specific Analysis
================================

.. image:: https://img.shields.io/pypi/v/wasp2
   :target: https://pypi.org/project/wasp2/
   :alt: PyPI

.. image:: https://github.com/Jaureguy760/WASP2-exp/workflows/WASP2%20Tests/badge.svg
   :target: https://github.com/Jaureguy760/WASP2-exp/actions
   :alt: Tests

WASP2 is a comprehensive suite of tools for unbiased allele-specific analysis of next-generation sequencing data. It addresses reference bias in read mapping and provides statistical methods for detecting allelic imbalance.

Features
--------

* **Unbiased Mapping**: WASP algorithm for correcting reference bias
* **Allele Counting**: Count allele-specific reads from BAM files
* **Statistical Analysis**: Beta-binomial models for allelic imbalance detection
* **Single-Cell Support**: Specialized tools for single-cell RNA-seq
* **Type-Safe**: 100% type hint coverage for robust code
* **Well-Tested**: Comprehensive regression and integration tests

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

   tutorials/scrna_seq
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
