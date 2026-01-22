I/O Module API
==============

The I/O module provides efficient variant data loading from multiple formats including VCF, BCF, and PGEN.

Variant Source Base
-------------------

.. automodule:: wasp2.io.variant_source
   :members:
   :undoc-members:
   :show-inheritance:

VCF Source
----------

Standard VCF/BCF file reading using pysam.

.. automodule:: wasp2.io.vcf_source
   :members:
   :undoc-members:
   :show-inheritance:

cyvcf2 Source
-------------

High-performance VCF reading using cyvcf2 (7x faster parsing).

.. automodule:: wasp2.io.cyvcf2_source
   :members:
   :undoc-members:
   :show-inheritance:

PGEN Source
-----------

PLINK 2.0 PGEN format support (up to 25x faster I/O).

.. automodule:: wasp2.io.pgen_source
   :members:
   :undoc-members:
   :show-inheritance:

Compatibility Utilities
-----------------------

.. automodule:: wasp2.io.compat
   :members:
   :undoc-members:
   :show-inheritance:

Module Initialization
---------------------

.. automodule:: wasp2.io
   :members:
   :undoc-members:
   :show-inheritance:
