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

Verification
------------

.. code-block:: bash

   wasp2-count --help
   wasp2-map --help
   wasp2-analyze --help
