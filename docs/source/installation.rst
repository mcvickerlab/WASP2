Installation
============

Via Bioconda (Recommended)
--------------------------

`Bioconda <https://bioconda.github.io>`_ installs WASP2 and **all** dependencies
(samtools, bcftools, bedtools, htslib) in one command. Requires
`miniforge <https://github.com/conda-forge/miniforge>`_.

.. code-block:: bash

   mamba install -c conda-forge -c bioconda wasp2

Or with conda:

.. code-block:: bash

   conda install -c conda-forge -c bioconda wasp2

Available for Linux (x86_64, aarch64) and macOS (Intel, Apple Silicon) with
Python 3.11-3.12.

Via PyPI
--------

.. code-block:: bash

   pip install wasp2

Pre-built wheels include the Rust extension and bundled htslib for Linux
(x86_64, aarch64) and macOS (Intel, Apple Silicon) with Python 3.10-3.13.

.. note::

   The PyPI package does not include samtools, bcftools, or bedtools.
   Install them separately before running WASP2:

   * On Ubuntu/Debian: ``sudo apt-get install bcftools bedtools samtools``
   * On macOS: ``brew install bcftools bedtools samtools``
   * Via conda: ``mamba install -c bioconda samtools bcftools bedtools``

Via Docker
----------

WASP2 is available as a multi-platform Docker image (linux/amd64 + linux/arm64)
with all dependencies pre-installed:

.. code-block:: bash

   docker pull ghcr.io/mcvickerlab/wasp2:1.4.0

   # Run a command
   docker run --rm -v /path/to/data:/data ghcr.io/mcvickerlab/wasp2:1.4.0 \
       wasp2-count count-variants /data/sample.bam /data/variants.vcf

   # Interactive shell
   docker run -it --rm -v /path/to/data:/data ghcr.io/mcvickerlab/wasp2:1.4.0 bash

The image includes samtools, bcftools, bedtools, and the Rust-accelerated backend.

Via Singularity/Apptainer
-------------------------

For HPC environments that don't support Docker:

.. code-block:: bash

   singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.4.0
   singularity exec wasp2.sif wasp2-count --help

Development Installation
------------------------

For contributing or building from source:

.. code-block:: bash

   git clone https://github.com/mcvickerlab/WASP2.git
   cd WASP2
   pixi install              # resolves all deps including Rust toolchain
   pixi run build            # builds the Rust extension
   pixi run test             # runs the test suite

`pixi <https://pixi.sh>`_ resolves Python, Rust toolchain, samtools, bcftools,
bedtools, and htslib automatically. No system packages required.

Compiling pgenlib
~~~~~~~~~~~~~~~~~

WASP2 optionally uses `pgenlib <https://github.com/chrchang/plink-ng/tree/master/2.0/Python>`_
for PLINK2 file I/O. This requires a C compiler:

* On Ubuntu/Debian: ``sudo apt-get install build-essential python3-dev``
* On macOS: ``xcode-select --install``
* On RHEL/Fedora: ``sudo dnf install gcc gcc-c++ python3-devel``

pgenlib is installed automatically via pip when you install WASP2.

Verification
------------

.. code-block:: bash

   wasp2-count --help
   wasp2-map --help
   wasp2-analyze --help
