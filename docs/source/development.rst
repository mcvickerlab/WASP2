Development Guide
=================

Contributing to WASP2
---------------------

We welcome contributions! This guide helps you get started.

Development Setup
-----------------

Clone Repository
~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/Jaureguy760/WASP2-exp
   cd WASP2-exp

Install Development Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install -e ".[dev]"

This installs:
* pytest (testing)
* mypy (type checking)
* black (code formatting)
* flake8 (linting)
* pre-commit (git hooks)

Install Pre-commit Hooks
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pre-commit install

Hooks run automatically before each commit:
* Black formatting
* Flake8 linting
* mypy type checking
* Quick tests

Code Standards
--------------

Type Hints
~~~~~~~~~~

WASP2 has 100% type hint coverage. All new code must include type hints:

.. code-block:: python

   def count_alleles(
       bam_file: str,
       vcf_file: str,
       min_count: int = 10
   ) -> pd.DataFrame:
       """Count alleles from BAM file."""
       ...

Formatting
~~~~~~~~~~

Use Black with 100-character lines:

.. code-block:: bash

   black src/ --line-length=100

Linting
~~~~~~~

Pass Flake8 checks:

.. code-block:: bash

   flake8 src/ --max-line-length=100

Testing
-------

Run Tests Locally
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # All tests
   pytest tests/ -v

   # Fast tests only (skip slow integration tests)
   pytest tests/ -v -m "not slow"

   # With coverage
   pytest tests/ --cov=src --cov-report=html

Test Requirements
~~~~~~~~~~~~~~~~~

* All new features need tests
* Maintain >80% code coverage
* Tests must pass in CI before merge

Type Checking
-------------

Run mypy:

.. code-block:: bash

   mypy src/counting/ src/mapping/ src/analysis/

All code must pass mypy with 0 errors.

CI/CD Pipeline
--------------

GitHub Actions
~~~~~~~~~~~~~~

Tests run automatically on every push:
* Python 3.10 and 3.11
* Type checking (mypy)
* Unit tests (pytest)
* Full pipeline validation
* Documentation build

CI must pass before PR can be merged.

Pre-commit Hooks
~~~~~~~~~~~~~~~~

Local checks before commit:
* Code formatting
* Type checking
* Quick tests

To bypass (not recommended):

.. code-block:: bash

   git commit --no-verify

Pull Request Process
--------------------

1. Fork & Branch
~~~~~~~~~~~~~~~~

.. code-block:: bash

   git checkout -b feature/my-feature

2. Develop & Test
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Make changes
   vim src/analysis/my_feature.py

   # Add type hints
   # Write tests
   # Run locally
   pytest tests/ -v
   mypy src/

3. Commit
~~~~~~~~~

.. code-block:: bash

   git add src/analysis/my_feature.py tests/test_my_feature.py
   git commit -m "Add my feature"

   # Pre-commit hooks run automatically

4. Push & PR
~~~~~~~~~~~~

.. code-block:: bash

   git push origin feature/my-feature

   # Open PR on GitHub
   # CI will run automatically
   # Request review

Code Review
-----------

PRs are reviewed for:
* Correctness
* Type safety
* Test coverage
* Documentation
* Code style

Project Structure
-----------------

.. code-block:: text

   WASP2-exp/
   ├── src/
   │   ├── counting/       # Allele counting
   │   ├── mapping/        # WASP remapping
   │   └── analysis/       # Statistical analysis
   ├── tests/
   │   └── regression/     # Regression tests
   ├── docs/               # Sphinx documentation
   ├── scripts/            # Utility scripts
   ├── baselines/          # Test baselines
   └── test_data/          # Example data

Building Documentation
----------------------

.. code-block:: bash

   cd docs
   make html
   open build/html/index.html

Documentation must build without warnings.

Release Process
---------------

1. Update version in ``pyproject.toml``
2. Update ``docs/source/changelog.rst``
3. Merge to main
4. Tag release: ``git tag v1.1.0``
5. Push tag: ``git push origin v1.1.0``
6. Publish to PyPI: ``python -m build && twine upload dist/*``

AI-Assisted Development
-----------------------

WASP2 pipeline development benefits from AI tooling. See the full integration guide:
:doc:`/seqera_ai_integration`

Recommended Workflow
~~~~~~~~~~~~~~~~~~~~

1. **Design**: Use Claude Code for architecture and complex logic
2. **Generate**: Use Seqera AI for DSL2 syntax and nf-test templates
3. **Validate**: Use Anthropic life-sciences scripts for environment checks
4. **Review**: Use Claude Code for code review and optimization

Tool Selection
~~~~~~~~~~~~~~

* **Architecture and design** → Claude Code
* **Nextflow DSL2 syntax** → Seqera AI
* **nf-test generation** → Seqera AI
* **Environment validation** → ``nextflow run . -profile test -preview``

Getting Help
------------

* **Issues**: https://github.com/Jaureguy760/WASP2-exp/issues
* **Discussions**: GitHub Discussions
* **Email**: Contact maintainers

License
-------

WASP2 is released under the MIT License. See LICENSE file.
