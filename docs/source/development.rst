Development Guide
=================

Contributions are welcome. This page covers the WASP2-specific bits; for
generic Python tooling (pytest, black, pre-commit) refer to each project's
own docs.

Setup
-----

.. code-block:: bash

   git clone https://github.com/mcvickerlab/WASP2
   cd WASP2
   pip install -e ".[dev]"
   pre-commit install

The ``[dev]`` extra pulls in pytest, mypy, black, flake8, pre-commit, and
Sphinx. Pre-commit runs formatting + linting + mypy + a quick test subset
on every commit.

Code standards
--------------

- **Type hints are required** on new public API. WASP2 maintains full
  mypy coverage on ``src/counting``, ``src/mapping``, and ``src/analysis``.
- **Line length 100**, enforced by black.
- **Ruff** runs as the primary lint tool; flake8 is retained for
  docstring-specific rules.
- **Tests** live in ``tests/``. All new features need a test. Maintain
  ≥80% line coverage.

.. code-block:: bash

   pytest tests/ -v                  # full suite
   pytest tests/ -v -m "not slow"    # fast subset
   pytest tests/ --cov=src           # coverage

   mypy src/counting/ src/mapping/ src/analysis/

Branching
---------

- ``master`` is the release branch. Tagged versions (e.g., ``v1.4.1``)
  trigger PyPI + Docker publish via GitHub Actions.
- ``dev`` is the integration branch. Feature PRs target ``dev``, not
  ``master``.
- Use short, kebab-case feature branches off ``dev``: ``feat/...``,
  ``fix/...``, ``docs/...``.
- After a batch of features merges to ``dev``, a ``dev → master`` PR
  promotes them; CI runs on both sides of that merge.

Pull requests
-------------

1. Branch off ``dev``.
2. Commit changes with descriptive messages; the squash-merge commit on
   ``dev`` is what ends up in the changelog, so the PR title and body
   should be clean.
3. CI must be green: ruff + Rust tests + pytest (Python 3.10 / 3.11 /
   3.12) + CodeQL.
4. Code review looks at correctness, type safety, test coverage, docs,
   and API surface.

Rust layer
----------

The ``wasp2_rust`` extension is built with maturin against Python 3.10+.

.. code-block:: bash

   cd rust/
   cargo fmt --check
   cargo clippy -- -D warnings
   cargo test
   maturin develop --release --manifest-path Cargo.toml

Changes to Rust code must keep the parity test
(``tests/test_rust_python_counting_parity.py``) green — the Rust counter
must produce byte-identical output to the Python reference implementation
on the shared chr21 fixture.

Project layout
--------------

.. code-block:: text

   WASP2/
   ├── src/
   │   ├── counting/     # allele counting (Python + Rust bridge)
   │   ├── mapping/      # WASP remap-and-filter
   │   └── analysis/     # beta-binomial LRT, FDR
   ├── rust/             # Rust extension (bam filter, counter)
   ├── tests/            # pytest + regression
   ├── docs/             # Sphinx source
   └── pipelines/        # Nextflow DSL2 workflows

Docs
----

.. code-block:: bash

   cd docs && make html

The live site is built and deployed by GitHub Pages on every push to
``master``. Docs must build without warnings; cross-references are
checked by ``sphinx-build -n`` in CI.

Releases
--------

1. Update ``pyproject.toml`` and ``rust/Cargo.toml`` versions (they share
   one source of truth — ``pyproject.toml`` reads from Cargo).
2. Update ``CHANGELOG.md``.
3. Merge ``dev → master``.
4. Tag: ``git tag vX.Y.Z && git push origin vX.Y.Z``.
5. The ``release.yml`` workflow builds wheels and publishes to PyPI via
   OIDC trusted publishing; ``docker.yml`` publishes the ghcr.io image.

Getting help
------------

- Issues and feature requests: https://github.com/mcvickerlab/WASP2/issues
- Discussions: GitHub Discussions
