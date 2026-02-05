# WASP2 v1.4.0 Release Audit Checklist

## Pre-Release Verification

Use this checklist to track audit progress. Each team member should complete their assigned section.

---

## Code Quality (Software Engineer)

### CI/CD & Build
- [ ] All CI checks pass on main branch
- [ ] CI passes on all platforms (Linux, macOS, Windows)
- [ ] Docker build completes successfully
- [ ] Singularity build completes successfully
- [ ] Wheel builds for all platforms (manylinux, macos, windows)

### Test Coverage
- [ ] Test coverage meets threshold (≥30%)
- [ ] Audit #200 test coverage gaps reviewed and prioritized
- [ ] No critical paths without test coverage
- [ ] Integration tests pass with real data

### Code Quality
- [ ] No security vulnerabilities (`bandit -r src/`)
- [ ] No dependency vulnerabilities (`pip-audit`)
- [ ] Version consistency verified (`scripts/check-version-consistency.sh`)
- [ ] Pre-commit hooks pass (`pre-commit run --all-files`)

### Rust Component
- [ ] `cargo test` passes all tests
- [ ] `cargo clippy` shows no warnings
- [ ] Rust-Python binding works correctly
- [ ] Memory safety verified (no unsafe blocks without justification)

**Sign-off:** __________________ Date: __________

---

## Scientific Accuracy (Bioinformatician)

### Statistical Model
- [ ] Beta-binomial model in `as_analysis.py` validated
- [ ] Dispersion clamping works correctly (Issue #228)
- [ ] P-value calculations verified
- [ ] Effect size calculations verified

### Benchmark Concordance
- [ ] GATK ASEReadCounter comparison (r² > 0.99)
- [ ] Sanity tests pass with chr21 data
- [ ] Edge cases tested (zero counts, single allele)

### Variant Handling
- [ ] SNP handling verified
- [ ] INDEL handling verified (`test_indel_correctness.py`)
- [ ] Multi-allelic variant handling verified
- [ ] VCF parsing robust

### Single-Cell Analysis
- [ ] Single-cell output format correct
- [ ] Cell barcode handling verified
- [ ] Aggregation logic validated

**Sign-off:** __________________ Date: __________

---

## Documentation (PI)

### Release Documentation
- [ ] CHANGELOG.md accurately reflects all changes
- [ ] README.md is up-to-date
- [ ] CITATION.cff has correct version and date
- [ ] Release notes drafted

### User Documentation
- [ ] Installation instructions work
- [ ] Tutorials are runnable
- [ ] API documentation complete
- [ ] Error messages documented

### Security & Compliance
- [ ] SECURITY_AUDIT.md reviewed and approved
- [ ] No sensitive data in repository
- [ ] License file present and correct
- [ ] Code of conduct present

### Final Approval
- [ ] All critical issues resolved
- [ ] All high severity issues resolved or deferred with justification
- [ ] Release candidate tested end-to-end

**Sign-off:** __________________ Date: __________

---

## Real-World Testing (Staff Research Scientist)

### Data Compatibility
- [ ] Tested with real CRISPR/ASE data
- [ ] VCF input handling verified with lab genotypes
- [ ] PGEN input handling verified (if applicable)
- [ ] BAM/CRAM handling verified

### Container Testing
- [ ] Docker container runs on lab systems
- [ ] Singularity container runs on HPC
- [ ] Memory usage acceptable
- [ ] Runtime acceptable

### Output Usability
- [ ] Output files are interpretable
- [ ] Output can be used for downstream analysis
- [ ] Visualization recommendations work
- [ ] Results match expectations from known samples

### Documentation Clarity
- [ ] Non-computational users can follow tutorials
- [ ] Error messages are actionable
- [ ] Help text is clear
- [ ] Examples are relevant to real use cases

**Sign-off:** __________________ Date: __________

---

## Issue Summary

### Critical Issues (Must Fix)
| Issue # | Description | Status | Assignee |
|---------|-------------|--------|----------|
| | | | |

### High Severity Issues (Should Fix)
| Issue # | Description | Status | Assignee |
|---------|-------------|--------|----------|
| | | | |

### Deferred Issues (Document for Future)
| Issue # | Description | Reason for Deferral |
|---------|-------------|---------------------|
| | | |

---

## Final Release Checklist

### Pre-Release
- [ ] All critical issues closed
- [ ] All high severity issues closed or documented
- [ ] Audit sign-offs complete from all team members
- [ ] Release candidate tag created
- [ ] Final test pass on release candidate

### Release
- [ ] Version number updated
- [ ] CHANGELOG finalized
- [ ] Git tag created and pushed
- [ ] PyPI release published
- [ ] Conda release published (bioconda PR)
- [ ] GitHub release created
- [ ] Documentation deployed

### Post-Release
- [ ] Announcement posted
- [ ] Known issues documented
- [ ] Next version roadmap updated

---

## Approval Signatures

| Role | Name | Signature | Date |
|------|------|-----------|------|
| PI | | | |
| Software Engineer | | | |
| Bioinformatician | | | |
| Staff Research Scientist | | | |

**Release Approved:** Yes / No

**Release Version:** v1.4.0

**Release Date:** __________
