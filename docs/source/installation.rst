Installation
============

Requirements
------------

System Dependencies
~~~~~~~~~~~~~~~~~~~

WASP2 requires:

* bcftools >= 1.10
* bedtools >= 2.29
* samtools >= 1.10

On Ubuntu/Debian:

.. code-block:: bash

   sudo apt-get install bcftools bedtools samtools

On macOS with Homebrew:

.. code-block:: bash

   brew install bcftools bedtools samtools

Python Requirements
~~~~~~~~~~~~~~~~~~~

* Python >= 3.10
* See pyproject.toml for full list

Installation
------------

Via PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install wasp2

Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/Jaureguy760/WASP2-exp
   cd WASP2-exp
   pip install -e ".[dev]"

Conda Installation
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   conda env create -f environment.yml
   conda activate wasp2

Docker Installation
~~~~~~~~~~~~~~~~~~~

WASP2 is available as a Docker image with all dependencies pre-installed.
This is the easiest way to get started, especially on systems where installing
bioinformatics tools is challenging.

**Pull from GitHub Container Registry:**

.. code-block:: bash

   docker pull ghcr.io/jaureguy760/wasp2:latest

**Run WASP2 commands:**

.. code-block:: bash

   # Run counting
   docker run -v /path/to/data:/data ghcr.io/jaureguy760/wasp2:latest \
       wasp2-count count-variants /data/sample.bam /data/variants.vcf

   # Interactive shell
   docker run -it -v /path/to/data:/data ghcr.io/jaureguy760/wasp2:latest bash

**Build locally (optional):**

.. code-block:: bash

   git clone https://github.com/Jaureguy760/WASP2-final
   cd WASP2-final
   docker build -t wasp2:local .

For detailed Docker usage including GPU support and Singularity conversion,
see the `Container Usage Guide <https://github.com/Jaureguy760/WASP2-final/blob/main/docs/CONTAINER_USAGE.md>`_.

Singularity/Apptainer
~~~~~~~~~~~~~~~~~~~~~

For HPC environments that don't support Docker, use Singularity/Apptainer:

.. code-block:: bash

   # Pull from GitHub Container Registry
   singularity pull wasp2.sif docker://ghcr.io/jaureguy760/wasp2:latest

   # Run WASP2 commands
   singularity exec wasp2.sif wasp2-count --help

   # Build from definition file
   singularity build wasp2.sif Singularity.def

Verification
------------

.. code-block:: bash

   wasp2-count --help
   wasp2-map --help
   wasp2-analyze --help
