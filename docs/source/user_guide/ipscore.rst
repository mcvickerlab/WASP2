iPSCORE Utilities
=================

Overview
--------

``wasp2-ipscore`` provides helper commands for working with the iPSCORE
multi-tissue allelic imbalance resource bundled with WASP2.

Commands
--------

``inventory``
~~~~~~~~~~~~~

Validate the expected iPSCORE data inventory and optionally write a TSV report.

.. code-block:: bash

   wasp2-ipscore inventory --output inventory.tsv

``manifest``
~~~~~~~~~~~~

Create a unified sample manifest across tissues and assays.

.. code-block:: bash

   wasp2-ipscore manifest --output manifest.csv --format csv

``qtls``
~~~~~~~~

Load and summarize QTL resources, optionally filtering by tissue.

.. code-block:: bash

   wasp2-ipscore qtls --tissue iPSC --output qtls.tsv

``validate``
~~~~~~~~~~~~

Run the combined inventory, manifest, and QTL validation workflow.

.. code-block:: bash

   wasp2-ipscore validate

Notes
-----

These commands are resource-management utilities. They do not replace the main
``wasp2-map``, ``wasp2-count``, or ``wasp2-analyze`` analysis workflows.
