Installation
============

Published Packages
------------------

Package registries and container tags may lag behind ``dev``. This source
branch targets the unreleased 1.5.0 version; installing a published package
does not necessarily include its changes. Check ``wasp2-count --version``
and use the source installation below when testing development fixes.

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

Select an existing published tag from the
`container registry <https://github.com/mcvickerlab/WASP2/pkgs/container/wasp2>`_.
Replace ``PUBLISHED_TAG`` below with that tag. A source version does not
guarantee that a matching container has been published.

.. code-block:: bash

   docker pull ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG

   docker run --rm ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG wasp2-count --help
   docker run --rm ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG wasp2-map --help
   docker run --rm ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG wasp2-analyze --help

Via Singularity/Apptainer
-------------------------

For HPC environments that require SIF images:

.. code-block:: bash

   singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG
   singularity exec wasp2.sif wasp2-count --help

or with Apptainer:

.. code-block:: bash

   apptainer pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG
   apptainer exec wasp2.sif wasp2-count --help

Development Installation
------------------------

Source builds require Rust 1.88 or newer, a C/C++ compiler, and libclang.
The build backend declares CMake as a build dependency. Pixi supplies these
tools together with the runtime dependencies; published wheels do not need
the compilation toolchain.

.. code-block:: bash

   git clone --branch dev https://github.com/mcvickerlab/WASP2.git
   cd WASP2
   pixi install
   pixi run verify

Verification
------------

.. code-block:: bash

   wasp2-count --help
   wasp2-map --help
   wasp2-analyze --help
