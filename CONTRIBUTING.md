# Contributing to WASP2

Thank you for your interest in contributing to WASP2! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [First-Time Contributors](#first-time-contributors)
- [Development Setup](#development-setup)
- [Running Tests](#running-tests)
- [Code Style](#code-style)
- [Branch Workflow](#branch-workflow)
- [Issue Guidelines](#issue-guidelines)
- [Pull Request Process](#pull-request-process)
- [Troubleshooting](#troubleshooting)
- [License](#license)

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## First-Time Contributors

New to WASP2? Welcome! Here's how to get started:

1. Look for issues labeled [`good first issue`](https://github.com/Jaureguy760/WASP2-final/labels/good%20first%20issue) or [`help wanted`](https://github.com/Jaureguy760/WASP2-final/labels/help%20wanted)
2. Comment on the issue to let maintainers know you're working on it
3. Fork the repository and follow the [Development Setup](#development-setup)
4. Submit your PR following the [Pull Request Process](#pull-request-process)

We welcome contributions of all sizes, from typo fixes to new features!

## Development Setup

WASP2 is a hybrid Python/Rust project. You'll need both toolchains to build from source.

### Prerequisites

| Requirement | Version | Notes |
|-------------|---------|-------|
| Python | 3.10+ | 3.10, 3.11, 3.12 supported (conda env uses 3.11) |
| Rust | 1.70+ | Edition 2021 |
| Git | 2.0+ | For version control |
| C compiler | gcc/clang | Required for native extensions |

> **Note:** WASP2 is developed and tested on Linux and macOS. Windows is not officially supported.

### Setting Up the Environment

1. **Clone the repository:**

   ```bash
   git clone https://github.com/Jaureguy760/WASP2-final.git
   cd WASP2-final
   ```

2. **Create a Python virtual environment:**

   Using conda (recommended):
   ```bash
   conda env create -f environment.yml
   conda activate WASP2
   ```

   Or using venv:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Linux/macOS
   pip install -e ".[dev]"    # Includes maturin and dev tools
   ```

3. **Install the Rust toolchain** (if not already installed):

   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source ~/.cargo/env
   ```

4. **Build the Rust extension:**

   ```bash
   # Development build (faster compile, debug symbols)
   maturin develop -m rust/Cargo.toml

   # Release build (optimized)
   maturin develop --release -m rust/Cargo.toml
   ```

5. **Install pre-commit hooks:**

   ```bash
   pip install pre-commit
   pre-commit install
   ```

6. **Verify the installation:**

   ```bash
   make verify-cli
   ```

### Makefile Commands

The project includes a Makefile with common development commands:

```bash
make help          # Show all available commands
make build         # Build Rust extension and install package
make rust-dev      # Build Rust extension in debug mode
make test          # Run all tests
make lint          # Run all linters
make format        # Format all code
make security      # Run security checks (bandit, cargo audit)
```

## Running Tests

### Python Tests

```bash
# Run all tests
pytest tests/ -v

# Run quick validation tests
pytest tests/test_validation_quick.py -v

# Run specific test markers
pytest tests/ -v -m "unit"           # Unit tests only
pytest tests/ -v -m "integration"    # Integration tests only
pytest tests/ -v -m "rust"           # Rust backend tests
pytest tests/ -v -m "not slow"       # Exclude slow tests

# Run with coverage
pytest tests/ --cov=src --cov-report=html
```

### Rust Tests

```bash
cd rust
cargo test
```

### Sanity Tests (Real Data)

For contributors with access to test data:

```bash
# Download test data from GitHub release
make download-sanity-data

# Run sanity tests
make test-sanity
```

## Code Style

### Python

We use **[ruff](https://github.com/astral-sh/ruff)** for linting and formatting:

```bash
# Check for issues
ruff check src/ tests/

# Auto-fix issues
ruff check --fix src/ tests/

# Format code
ruff format src/ tests/
```

Key style settings (configured in `pyproject.toml`):
- Line length: 100 characters
- Target Python version: 3.10
- Import sorting: isort-compatible via ruff

### Rust

We use **cargo fmt** and **clippy**:

```bash
cd rust

# Format code
cargo fmt

# Run linter
cargo clippy -- -D warnings
```

### Pre-commit Hooks

Pre-commit hooks run automatically on `git commit`. To run manually:

```bash
pre-commit run --all-files
```

The hooks include:
- **ruff**: Python linting and formatting
- **pre-commit-hooks**: File hygiene (trailing whitespace, end-of-file, YAML validation, large files, merge conflicts, private keys, AST validation)
- **bandit**: Python security linting
- **gitleaks**: Secret detection
- **basedpyright**: Type checking

### Security Scanning

Run security checks before submitting PRs:

```bash
# Python security audit
bandit -c pyproject.toml -r src/

# Dependency vulnerability scan
pip-audit

# Rust security audit
cd rust && cargo audit
```

## Branch Workflow

We use a feature branch workflow with PRs to `main`:

```
feature/* ─┐
fix/*     ─┼──→ main
docs/*    ─┘
```

### Branch Naming

- `feature/<description>` - New features
- `fix/<description>` - Bug fixes
- `docs/<description>` - Documentation updates
- `refactor/<description>` - Code refactoring
- `test/<description>` - Test additions or fixes

### Workflow

1. **Create a feature branch from `main`:**

   ```bash
   git checkout main
   git pull origin main
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes and commit:**

   ```bash
   git add .
   git commit -m "feat: add your feature description"
   ```

3. **Push and create a pull request:**

   ```bash
   git push -u origin feature/your-feature-name
   ```

### Commit Messages

We follow [Conventional Commits](https://www.conventionalcommits.org/):

| Prefix | Purpose |
|--------|---------|
| `feat:` | New feature |
| `fix:` | Bug fix |
| `docs:` | Documentation changes |
| `style:` | Code style changes (formatting, no logic change) |
| `refactor:` | Code refactoring |
| `test:` | Adding or updating tests |
| `chore:` | Maintenance tasks |
| `perf:` | Performance improvements |
| `ci:` | CI/CD changes |

## Issue Guidelines

### Before Creating an Issue

1. Search existing issues to avoid duplicates
2. Check the [documentation](https://jaureguy760.github.io/WASP2-final/) for answers
3. Ensure you're using the latest version

### Bug Reports

Include:
- WASP2 version (`pip show wasp2` or `python -c "import wasp2; print(wasp2.__version__)"`)
- Python version (`python --version`)
- Operating system and version
- Minimal reproducible example
- Expected vs. actual behavior
- Full error traceback

### Feature Requests

Include:
- Clear description of the proposed feature
- Use case and motivation
- Example of how it would be used

## Pull Request Process

No Contributor License Agreement (CLA) is required. By submitting a PR, you agree your contributions will be licensed under the MIT License.

1. **Ensure your code passes all checks:**

   ```bash
   make lint
   make test
   ```

2. **Update documentation** if needed (docstrings, README, etc.)

3. **Add tests** for new functionality

4. **Create the pull request:**
   - Use a clear, descriptive title
   - Reference any related issues (e.g., "Fixes #123")
   - Describe what changes were made and why
   - Include any relevant screenshots or output

5. **Address review feedback** promptly

### PR Checklist

- [ ] Code follows the project's style guidelines
- [ ] Tests pass locally (`make test`)
- [ ] Linting passes (`make lint`)
- [ ] Security checks pass (`make security`)
- [ ] New code is tested
- [ ] Documentation is updated (if applicable)
- [ ] Commit messages follow conventional commits

## Troubleshooting

### Common Issues

**Maturin build fails with "cargo not found":**
```bash
source ~/.cargo/env  # Add Rust to PATH
```

**Pre-commit hooks fail on first run:**
```bash
pre-commit run --all-files  # Run once to cache hooks
```

**Import errors after building:**
```bash
pip install -e ".[dev]" --force-reinstall
```

**Rust extension not loading:**
```bash
maturin develop --release -m rust/Cargo.toml
python -c "import wasp2_rust; print('OK')"
```

## License

By contributing to WASP2, you agree that your contributions will be licensed under the [MIT License](LICENSE).

---

Thank you for contributing to WASP2!
