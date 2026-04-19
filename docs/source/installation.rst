Installation
============

Via Bioconda
------------

.. code-block:: bash

   mamba install -c conda-forge -c bioconda wasp2

Bioconda installs WASP2 together with the external command-line dependencies
required by the workflows, including ``samtools``, ``bcftools``, and
``bedtools``.

Via PyPI
--------

.. code-block:: bash

   pip install wasp2

The PyPI wheel includes the WASP2 Python package and Rust extension, but it
does not install external tools such as ``samtools``, ``bcftools``, or
``bedtools``. Install those separately before running mapping or counting.

Via Docker
----------

The Docker image is the most reproducible fully bundled option available in
this environment.

.. code-block:: bash

   docker pull ghcr.io/mcvickerlab/wasp2:1.4.1

   docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.1 wasp2-count --help
   docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.1 wasp2-map --help
   docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.1 wasp2-analyze --help

Via Singularity/Apptainer
-------------------------

For HPC environments that require SIF images:

.. code-block:: bash

   singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.4.1
   singularity exec wasp2.sif wasp2-count --help

or with Apptainer:

.. code-block:: bash

   apptainer pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.4.1
   apptainer exec wasp2.sif wasp2-count --help

.. note::

   Docker was validated in this development environment. ``apptainer`` /
   ``singularity`` binaries were not available locally during this doc update,
   so those examples reflect the intended pull/exec workflow but were not
   executed here.

Development Installation
------------------------

.. code-block:: bash

   git clone https://github.com/mcvickerlab/WASP2.git
   cd WASP2
   pixi install
   pixi run verify

Verification
------------

.. code-block:: bash

   wasp2-count --help
   wasp2-map --help
   wasp2-analyze --help
