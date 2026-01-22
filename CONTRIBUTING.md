# Contributing to WASP2

Thank you for your interest in contributing to WASP2! This document provides guidelines for contributing.

## Development Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/Jaureguy760/WASP2-exp.git
   cd WASP2-exp
   ```

2. **Create conda environment**
   ```bash
   conda env create -f environment.yml
   conda activate WASP2
   ```

3. **Build the Rust extension**
   ```bash
   export LIBCLANG_PATH=$CONDA_PREFIX/lib
   export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
   export BINDGEN_EXTRA_CLANG_ARGS="-I/usr/include"
   maturin develop --release -m rust/Cargo.toml
   ```

4. **Install pre-commit hooks**
   ```bash
   pip install pre-commit
   pre-commit install
   ```

## Code Style

- **Python**: We use `black` for formatting and `flake8` for linting
- **Rust**: Use `cargo fmt` and `cargo clippy`
- Run `pre-commit run --all-files` before committing

## Testing

Run the test suite:
```bash
pytest tests/
```

Run validation against baselines:
```bash
export PYTHONPATH=$PWD
python validation/generate_baselines.py
python validation/compare_to_baseline.py
```

## Pull Request Process

1. Fork the repository and create a feature branch
2. Make your changes with clear, descriptive commits
3. Ensure all tests pass and pre-commit hooks succeed
4. Update documentation if needed
5. Submit a PR with a clear description of changes

## Reporting Issues

When reporting bugs, please include:
- WASP2 version (`pip show wasp2`)
- Python version
- Operating system
- Minimal reproducible example
- Full error traceback

## Code of Conduct

Be respectful and constructive in all interactions. We're building software to help researchers - let's keep it collaborative!

## Questions?

Open an issue or reach out to the maintainers.
