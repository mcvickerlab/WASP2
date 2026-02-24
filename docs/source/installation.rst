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

Compiling pgenlib
~~~~~~~~~~~~~~~~~

WASP2 uses `pgenlib <https://github.com/chrchang/plink-ng/tree/master/2.0/Python>`_ for
efficient PLINK2 file I/O. This library requires compilation from source and needs
a C compiler.

**Prerequisites:**

On Ubuntu/Debian:

.. code-block:: bash

   sudo apt-get install build-essential python3-dev

On macOS:

.. code-block:: bash

   xcode-select --install  # Installs Command Line Tools with clang

On RHEL/CentOS/Fedora:

.. code-block:: bash

   sudo dnf install gcc gcc-c++ python3-devel

**Installation:**

pgenlib is installed automatically via pip when you install WASP2:

.. code-block:: bash

   pip install pgenlib>=0.90

**Troubleshooting:**

If you encounter compilation errors:

1. **Missing Python headers**: Install ``python3-dev`` (Debian/Ubuntu) or ``python3-devel`` (RHEL/Fedora)
2. **No C compiler**: Install ``build-essential`` (Debian/Ubuntu) or ``gcc`` (RHEL/Fedora)
3. **macOS errors**: Ensure Xcode Command Line Tools are installed: ``xcode-select --install``
4. **Conda environments**: The environment.yml already includes Rust and Clang for PyO3 compilation

If compilation still fails, use the Docker image which has pgenlib pre-installed:

.. code-block:: bash

   docker pull ghcr.io/mcvickerlab/wasp2:latest

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

   git clone https://github.com/mcvickerlab/WASP2.git
   cd WASP2
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

   docker pull ghcr.io/mcvickerlab/wasp2:latest

**Run WASP2 commands:**

.. code-block:: bash

   # Run counting
   docker run -v /path/to/data:/data ghcr.io/mcvickerlab/wasp2:latest \
       wasp2-count count-variants /data/sample.bam /data/variants.vcf

   # Interactive shell
   docker run -it -v /path/to/data:/data ghcr.io/mcvickerlab/wasp2:latest bash

**Build locally (optional):**

.. code-block:: bash

   git clone https://github.com/mcvickerlab/WASP2
   cd WASP2
   docker build -t wasp2:local .

For detailed Docker usage including GPU support and Singularity conversion,
see the `Container Usage Guide <https://github.com/mcvickerlab/WASP2/blob/main/docs/CONTAINER_USAGE.md>`_.

Singularity/Apptainer
~~~~~~~~~~~~~~~~~~~~~

For HPC environments that don't support Docker, use Singularity/Apptainer:

.. code-block:: bash

   # Pull from GitHub Container Registry
   singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:latest

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
