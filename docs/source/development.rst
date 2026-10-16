Development Guide
=====================

Setup
---------

WASP2 uses Python 3.10–3.13 and Rust 1.88 or newer. Install a C compiler and
the external tools required by the workflows: samtools, bcftools, and bedtools.

.. code-block:: bash

   git clone --branch dev https://github.com/mcvickerlab/WASP2.git
   cd WASP2
   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev,docs]"
   pre-commit install

The ``dev`` extra provides pytest, Ruff, basedpyright, mypy, build, and
security tools. Documentation dependencies are in the separate ``docs`` extra.
Alternatively, use the locked Pixi environment with ``pixi install`` and
``pixi run verify``.

Checks
----------

Ruff handles Python linting and formatting with a 100-character line limit.
Pre-commit also runs file checks, Bandit, Gitleaks, and basedpyright. Add type
annotations to new public APIs and meaningful regressions for behavior changes.

.. code-block:: bash

   ruff check src tests scripts
   ruff format --check src tests scripts
   basedpyright
   pytest tests/ -v
   pytest tests/ -v -m "not slow"
   cargo fmt --manifest-path rust/Cargo.toml --check
   cargo clippy --manifest-path rust/Cargo.toml -- -D warnings
   cargo test --manifest-path rust/Cargo.toml

Rebuild the extension after Rust changes:

.. code-block:: bash

   maturin develop --release --manifest-path rust/Cargo.toml

Parity tests compare native counting and bulk inference against explicit
Python references. Preserve the statistical contract and keep legacy precision
limitations distinct from corrected wide-integer reference results. See the
existing counting and analysis tests for fixture requirements and tolerances.

Project Layout
------------------

* ``src/counting`` — allele counting and the Rust bridge.
* ``src/mapping`` — WASP remap-and-filter commands.
* ``src/analysis`` — bulk Rust dispatch and Python/SciPy single-cell inference.
* ``rust`` — native extension and numerical routines.
* ``tests`` — unit, regression, and integration tests.
* ``docs`` — Sphinx documentation.
* ``pipelines`` — Nextflow workflows and their documentation.

Documentation
-----------------

.. code-block:: bash

   make -C docs html

The documentation workflow builds pushes and pull requests to ``dev`` and
``master``. GitHub Pages deployment uses ``master``. Update examples alongside
public interface changes and check generated pages and cross-references.

Pull Requests and Releases
------------------------------

Branch from ``dev`` and target ``dev`` for feature and fix pull requests.
Use descriptive titles and explain behavior changes and validation. Required
checks include Python 3.10–3.13, Rust, lint, documentation, and security checks.

Version metadata comes from ``rust/Cargo.toml`` via maturin. Before a release,
run ``python scripts/release_version.py --check``, align the changelog and
release metadata, and verify package/container build results. Promotion to
``master`` and publication are separate release operations; an updated source
version does not establish that registry packages or images exist.

See `CONTRIBUTING.md <https://github.com/mcvickerlab/WASP2/blob/dev/CONTRIBUTING.md>`_
for contribution details. Report problems using the
`issue forms <https://github.com/mcvickerlab/WASP2/issues/new/choose>`_.
