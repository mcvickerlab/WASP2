WASP2: Allele-Specific Analysis
================================

.. image:: https://img.shields.io/pypi/v/wasp2
   :target: https://pypi.org/project/wasp2/
   :alt: PyPI

.. image:: https://github.com/Jaureguy760/WASP2-exp/workflows/WASP2%20Tests/badge.svg
   :target: https://github.com/Jaureguy760/WASP2-exp/actions
   :alt: Tests

.. image:: https://img.shields.io/badge/python-3.10+-blue
   :alt: Python 3.10+

.. image:: https://img.shields.io/badge/rust-1.70+-orange
   :alt: Rust 1.70+

.. image:: https://img.shields.io/badge/license-MIT-green
   :target: https://github.com/Jaureguy760/WASP2-exp/blob/master/LICENSE
   :alt: License

WASP2 is a comprehensive suite of tools for unbiased allele-specific analysis of next-generation sequencing data. It addresses reference bias in read mapping and provides statistical methods for detecting allelic imbalance.

Quick Start
-----------

Install via pip:

.. code-block:: bash

   pip install wasp2

Count alleles from a BAM file:

.. code-block:: bash

   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples NA12878 \
       --out_file counts.tsv

Analyze allelic imbalance:

.. code-block:: bash

   wasp2-analyze find-imbalance counts.tsv --out_file results.tsv

Key Features
------------

**Allele-Specific Analysis**
    * **Count Module**: Quantify ref/alt allele reads at heterozygous SNPs
    * **Analysis Module**: Beta-binomial statistical testing for allelic imbalance
    * **Mapping Module**: WASP algorithm for unbiased read mapping

**Performance**
    * **Rust Acceleration**: Core algorithms implemented in Rust (10-25x faster)
    * **Multi-Format Support**: VCF, BCF, PGEN (up to 25x faster I/O)
    * **High-Performance VCF**: Optional cyvcf2 backend (7x faster parsing)

**Applications**
    * RNA-seq allele-specific expression (ASE)
    * ATAC-seq allelic chromatin accessibility
    * Single-cell RNA-seq/ATAC-seq
    * ChIP-seq allelic binding analysis

**Quality**
    * 100% type hint coverage for robust code
    * Comprehensive regression and integration tests
    * Well-documented API with examples

Documentation
-------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index
   tutorials/concepts
   tutorials/basic_workflow
   tutorials/rnaseq_ase
   tutorials/atacseq_ase
   tutorials/single_cell

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/counting
   user_guide/mapping
   user_guide/analysis

.. toctree::
   :maxdepth: 2
   :caption: CLI Reference

   cli/index
   cli/wasp2_count
   cli/wasp2_map
   cli/wasp2_analyze

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/counting
   api/mapping
   api/analysis
   api/io

.. toctree::
   :maxdepth: 2
   :caption: Data Formats

   data_formats/index
   data_formats/input_formats
   data_formats/output_formats

.. toctree::
   :maxdepth: 1
   :caption: Reference

   faq
   troubleshooting
   development
   changelog

Comparison with Other Tools
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 18 18 18 18

   * - Feature
     - WASP2
     - GATK ASEReadCounter
     - phASER
     - MBASED
   * - Mapping Bias Correction
     - Yes (WASP)
     - No
     - No
     - No
   * - Statistical Testing
     - Beta-binomial
     - No
     - Phasing only
     - Beta-binomial
   * - Single-Cell Support
     - Yes
     - No
     - No
     - No
   * - Performance
     - Fast (Rust)
     - Slow
     - Medium
     - Medium
   * - Variant Formats
     - VCF/BCF/PGEN
     - VCF only
     - VCF only
     - VCF only
   * - Indel Support
     - Yes
     - Yes
     - No
     - No

Citation
--------

If you use WASP2 in your research, please cite:

.. code-block:: bibtex

    @article{wasp2_2025,
      title={WASP2: High-performance allele-specific analysis of next-generation sequencing data},
      author={Ho, Aaron and Jaureguy, Jeff and McVicker, Graham},
      journal={Bioinformatics},
      year={2025}
    }

**Original WASP paper:**

van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015). WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nature Methods* 12:1061-1063. `doi:10.1038/nmeth.3582 <https://doi.org/10.1038/nmeth.3582>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
