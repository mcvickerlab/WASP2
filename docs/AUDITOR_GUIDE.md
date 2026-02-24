# WASP2 Release Audit Guide

## For Human Reviewers of AI-Generated Code

This guide helps audit team members review WASP2 code for the v1.4.0 release. WASP2 contains AI-generated code that requires human validation for scientific accuracy, code quality, and security.

---

## Getting Started

### 1. Accept GitHub Invitation
Check your email for a collaborator invitation to `mcvickerlab/WASP2`.

### 2. Clone the Repository
```bash
git clone https://github.com/mcvickerlab/WASP2.git
cd WASP2
```

### 3. Set Up Your Environment
```bash
# Option A: Using conda (recommended)
conda env create -f environment.yml
conda activate wasp2

# Option B: Using pip
pip install -e ".[dev]"
```

### 4. Verify Installation
```bash
# Run quick sanity check
pytest tests/sanity/ -v --tb=short

# Check the CLI works
wasp2 --help
```

---

## How to Report Findings

### Creating an Issue

1. Go to **Issues** → **New Issue**
2. Select the appropriate template:
   - **Audit Finding**: General code quality, security, documentation issues
   - **Scientific Validation**: Domain-specific scientific concerns
3. Fill in all required fields
4. Add appropriate labels if needed

### Issue Templates Quick Reference

| Template | Use When | Who Should Use |
|----------|----------|----------------|
| Audit Finding | Code bugs, security issues, test gaps | All team members |
| Scientific Validation | Statistical accuracy, biological correctness | Bioinformatician, Scientist |

---

## Severity Guidelines

| Severity | Definition | Examples |
|----------|------------|----------|
| **Critical** | Blocks release, data corruption risk | Statistical model produces wrong p-values; security vulnerability |
| **High** | Must fix before release, significant impact | Missing edge case that affects 10%+ of use cases |
| **Medium** | Should fix, moderate impact | Poor error messages; documentation unclear |
| **Low** | Nice to have, minor impact | Code style improvements; minor documentation typos |

---

## Your Review Focus (by Role)

### PI (Admin)
- Final release authorization
- CHANGELOG.md completeness and accuracy
- CITATION.cff verification
- Documentation quality review
- Security audit sign-off

**Key Files:**
- `CHANGELOG.md`
- `CITATION.cff`
- `SECURITY_AUDIT.md`
- `README.md`
- `docs/source/installation.rst`

### Software Engineer (Maintain)
- CI/CD validation
- Test coverage analysis
- Rust/Python integration correctness
- Build system and packaging
- Performance optimization

**Key Files:**
- `.github/workflows/`
- `tests/` (coverage gaps from audit #200)
- `rust/` (Rust implementation)
- `pyproject.toml`
- `Dockerfile`, `Singularity.def`

**Commands to Run:**
```bash
# Full test suite
pytest tests/ -v --tb=short

# Test coverage report
pytest tests/ --cov=wasp2 --cov-report=html

# Rust tests
cd rust && cargo test

# Version consistency
scripts/check-version-consistency.sh

# Security scan
pip-audit
bandit -r src/
```

### Bioinformatician (Write)
- Statistical model validation
- Benchmark concordance (vs GATK)
- Output format correctness
- Scientific documentation accuracy

**Key Files:**
- `src/analysis/as_analysis.py` (beta-binomial model)
- `src/counting/count_alleles.py` (allele counting logic)
- `tests/sanity/` (validation tests)
- `tests/test_indel_correctness.py`
- `benchmarking/`

**Validation Steps:**
```bash
# Run sanity tests with chr21 data
pytest tests/sanity/ -v

# Compare with GATK ASEReadCounter
# Expected: r² > 0.99
python benchmarking/compare_gatk.py

# Check dispersion clamping (Issue #228)
python -c "from wasp2.analysis import as_analysis; print(as_analysis.DISPERSION_BOUNDS)"
```

### Staff Research Scientist (Write)
- Real-world data testing
- VCF/PGEN input handling
- Output interpretability
- Documentation clarity for non-computational users

**Key Files:**
- `docs/source/tutorials/`
- `tests/data/` (test data validation)
- `src/io/` (input handling)
- `tutorials/`

**Testing Focus:**
- Run with real CRISPR/ASE data from lab
- Verify VCF handling with lab-generated genotypes
- Check Docker container on lab computing resources
- Validate outputs are usable for wet lab follow-up

---

## Key Files to Review

| Area | Files | What to Check |
|------|-------|---------------|
| **Statistics** | `src/analysis/as_analysis.py` | Beta-binomial model correctness |
| **Counting** | `src/counting/count_alleles.py` | Allele counting accuracy |
| **Tests** | `tests/` | Coverage gaps (see audit #200) |
| **Rust Core** | `rust/src/` | Memory safety, performance |
| **CI/CD** | `.github/workflows/` | All platforms pass |
| **Docs** | `docs/source/` | Clarity for new users |
| **Security** | `SECURITY_AUDIT.md` | Known issues addressed |

---

## Common Issues to Look For

### Scientific Accuracy
- [ ] Statistical assumptions match published methods
- [ ] Edge cases handled (zero counts, missing data)
- [ ] Dispersion parameters within valid bounds
- [ ] Output values in expected ranges

### Code Quality
- [ ] Error messages are informative
- [ ] Input validation catches invalid data
- [ ] No hardcoded paths or magic numbers
- [ ] Logging provides useful debugging info

### Security
- [ ] No sensitive data in test fixtures
- [ ] Input sanitization for file paths
- [ ] Dependencies have no known vulnerabilities

### Documentation
- [ ] Installation instructions work
- [ ] Examples are accurate and runnable
- [ ] Parameters are documented
- [ ] Edge cases explained

---

## Communication

### Questions During Review
- Create a Discussion thread for questions
- Tag relevant team members using @username
- Use the issue templates for actual bugs

### Progress Updates
- Update your assigned issues as you work
- Mark issues as "in progress" when investigating
- Close issues when resolved with a summary

---

## Timeline

| Milestone | Date | Description |
|-----------|------|-------------|
| Audit Start | TBD | All team members have access |
| Initial Review | TBD +1 week | First pass complete, critical issues identified |
| Fix Period | TBD +2 weeks | Address critical and high severity issues |
| Final Review | TBD +3 weeks | Verify fixes, sign-off |
| Release | TBD +4 weeks | v1.4.0 published |

---

## Resources

- [WASP2 Documentation](https://wasp2.readthedocs.io/)
- [Original WASP Paper](https://doi.org/10.1038/nmeth.3582)
- [GitHub Issue Templates Guide](https://docs.github.com/en/communities/using-templates-to-encourage-useful-issues-and-pull-requests)
- [SECURITY_AUDIT.md](../SECURITY_AUDIT.md) - Known security considerations
- [CHANGELOG.md](../CHANGELOG.md) - Version history
