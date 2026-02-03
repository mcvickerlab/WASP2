# WASP2 Makefile
# Common development targets for building, testing, and benchmarking

.PHONY: all build install test test-quick test-sanity lint format clean help
.PHONY: download-sanity-data sanity-data-local rust-build rust-test

# Configuration
PYTHON ?= python
MATURIN ?= maturin
PYTEST ?= pytest
RUFF ?= ruff
CARGO ?= cargo

# Project paths
RUST_DIR := rust
SRC_DIR := src
TESTS_DIR := tests
SANITY_DATA_DIR := tests/sanity/data

# Sanity test data configuration
SANITY_VERSION := v1
SANITY_TARBALL := wasp2-sanity-chr21-$(SANITY_VERSION).tar.xz
SANITY_RELEASE_URL := https://github.com/Jaureguy760/WASP2-final/releases/download/v1.3.0/$(SANITY_TARBALL)

# Local sanity data path (for development)
LOCAL_SANITY_DATA := /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/sanity_test

# =============================================================================
# Main targets
# =============================================================================

all: build  ## Build the project (default)

build: rust-build install  ## Build Rust extension and install package

install:  ## Install package in development mode
	$(PYTHON) -m pip install -e ".[dev]" --no-build-isolation -q

# =============================================================================
# Rust build targets
# =============================================================================

rust-build:  ## Build Rust extension with maturin
	$(MATURIN) build --release -m $(RUST_DIR)/Cargo.toml
	$(PYTHON) -m pip install $(RUST_DIR)/target/wheels/*.whl --force-reinstall -q

rust-dev:  ## Build Rust extension in debug mode (faster compile)
	$(MATURIN) develop -m $(RUST_DIR)/Cargo.toml

rust-test:  ## Run Rust unit tests
	cd $(RUST_DIR) && $(CARGO) test

rust-bench:  ## Run Rust benchmarks
	cd $(RUST_DIR) && $(CARGO) bench

# =============================================================================
# Testing targets
# =============================================================================

test:  ## Run all tests (excluding benchmarks)
	$(PYTEST) $(TESTS_DIR) -v --tb=short \
		--ignore=$(TESTS_DIR)/benchmarks \
		-m "not benchmark"

test-quick:  ## Run quick validation tests only
	$(PYTEST) $(TESTS_DIR)/test_validation_quick.py -v --tb=short

test-rust:  ## Run Rust-specific tests
	$(PYTEST) $(TESTS_DIR) -v --tb=short -m "rust"

test-integration:  ## Run integration tests
	$(PYTEST) $(TESTS_DIR) -v --tb=short -m "integration"

test-sanity:  ## Run sanity tests with real chr21 data
	$(PYTEST) $(TESTS_DIR)/sanity -v --tb=short -x

test-all:  ## Run all tests including sanity and slow tests
	$(PYTEST) $(TESTS_DIR) -v --tb=short \
		--ignore=$(TESTS_DIR)/benchmarks

# =============================================================================
# Sanity data management
# =============================================================================

download-sanity-data:  ## Download sanity test data from GitHub release
	@echo "Downloading sanity data from $(SANITY_RELEASE_URL)..."
	@mkdir -p $(SANITY_DATA_DIR)
	@if command -v wget > /dev/null; then \
		wget -q -O $(SANITY_DATA_DIR)/$(SANITY_TARBALL) $(SANITY_RELEASE_URL); \
	else \
		curl -sL -o $(SANITY_DATA_DIR)/$(SANITY_TARBALL) $(SANITY_RELEASE_URL); \
	fi
	@echo "Extracting..."
	@cd $(SANITY_DATA_DIR) && tar -xJf $(SANITY_TARBALL) --strip-components=1
	@rm -f $(SANITY_DATA_DIR)/$(SANITY_TARBALL)
	@echo "Sanity data ready in $(SANITY_DATA_DIR)/"

sanity-data-local:  ## Link sanity data from local HPC path (development)
	@if [ -d "$(LOCAL_SANITY_DATA)" ]; then \
		mkdir -p $(SANITY_DATA_DIR); \
		ln -sf $(LOCAL_SANITY_DATA)/chr21.bam $(SANITY_DATA_DIR)/chr21.bam; \
		ln -sf $(LOCAL_SANITY_DATA)/chr21.bam.bai $(SANITY_DATA_DIR)/chr21.bam.bai; \
		ln -sf $(LOCAL_SANITY_DATA)/chr21.vcf.gz $(SANITY_DATA_DIR)/chr21.vcf.gz; \
		ln -sf $(LOCAL_SANITY_DATA)/chr21.vcf.gz.tbi $(SANITY_DATA_DIR)/chr21.vcf.gz.tbi; \
		ln -sf $(LOCAL_SANITY_DATA)/expected_counts.tsv $(SANITY_DATA_DIR)/expected_counts.tsv; \
		ln -sf $(LOCAL_SANITY_DATA)/expected_r1.fq.gz $(SANITY_DATA_DIR)/expected_r1.fq.gz; \
		ln -sf $(LOCAL_SANITY_DATA)/expected_r2.fq.gz $(SANITY_DATA_DIR)/expected_r2.fq.gz; \
		ln -sf $(LOCAL_SANITY_DATA)/expected_analysis.tsv $(SANITY_DATA_DIR)/expected_analysis.tsv; \
		echo "Linked sanity data from $(LOCAL_SANITY_DATA)"; \
	else \
		echo "Local sanity data not found at $(LOCAL_SANITY_DATA)"; \
		echo "Run 'make download-sanity-data' to download from GitHub"; \
		exit 1; \
	fi

clean-sanity-data:  ## Remove downloaded sanity test data
	rm -rf $(SANITY_DATA_DIR)

# =============================================================================
# Benchmarking targets
# =============================================================================

benchmark:  ## Run all benchmarks
	$(PYTHON) $(TESTS_DIR)/benchmarks/run_benchmarks.py

benchmark-quick:  ## Run quick benchmark subset
	$(PYTHON) $(TESTS_DIR)/benchmarks/run_benchmarks.py --quick

benchmark-figures:  ## Regenerate benchmark figures
	$(PYTHON) $(TESTS_DIR)/benchmarks/run_benchmarks.py --figures-only

benchmark-list:  ## List available benchmark groups
	$(PYTHON) $(TESTS_DIR)/benchmarks/run_benchmarks.py --list-groups

# =============================================================================
# Code quality targets
# =============================================================================

lint:  ## Run all linters
	$(RUFF) check $(SRC_DIR) $(TESTS_DIR)
	cd $(RUST_DIR) && $(CARGO) clippy -- -D warnings

format:  ## Format code
	$(RUFF) format $(SRC_DIR) $(TESTS_DIR)
	cd $(RUST_DIR) && $(CARGO) fmt

format-check:  ## Check code formatting without changes
	$(RUFF) format --check $(SRC_DIR) $(TESTS_DIR)
	cd $(RUST_DIR) && $(CARGO) fmt --check

typecheck:  ## Run type checking
	mypy $(SRC_DIR) --ignore-missing-imports

security:  ## Run security checks
	bandit -c pyproject.toml -r $(SRC_DIR)
	cd $(RUST_DIR) && cargo audit

# =============================================================================
# CLI verification
# =============================================================================

verify-cli:  ## Verify CLI tools are working
	wasp2-count --help > /dev/null && echo "wasp2-count: OK"
	wasp2-map --help > /dev/null && echo "wasp2-map: OK"
	wasp2-analyze --help > /dev/null && echo "wasp2-analyze: OK"
	$(PYTHON) -c "import wasp2_rust; print('wasp2_rust: OK')"

# =============================================================================
# Cleanup targets
# =============================================================================

clean:  ## Clean build artifacts
	rm -rf build/ dist/ *.egg-info
	rm -rf $(RUST_DIR)/target/wheels
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete

clean-all: clean clean-sanity-data  ## Clean everything including sanity data
	rm -rf $(RUST_DIR)/target
	rm -rf .pytest_cache .mypy_cache .ruff_cache

# =============================================================================
# Help
# =============================================================================

help:  ## Show this help message
	@echo "WASP2 Development Makefile"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
