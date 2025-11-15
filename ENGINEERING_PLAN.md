# WASP2-exp Engineering Plan
## Phase 1 & 2: Explore, Document, Clean & Refactor

**Goal**: Establish a solid, well-tested, well-documented foundation before adding new features (PGEN support, etc.)

---

## 📋 PHASE 1: EXPLORE & UNDERSTAND

### **1.1 Initial Code Inventory & Architecture Mapping**
**Duration**: 1-2 days
**Dependencies**: None
**Priority**: CRITICAL - Foundation for all other work

#### Tasks:
- [ ] **1.1.1** Map complete directory structure
  - Document purpose of each directory
  - Identify all entry points (bin/WASP2, module __main__.py files)
  - List all Python files with brief descriptions

- [ ] **1.1.2** Analyze external dependencies
  - Review environment.yml thoroughly
  - Document version constraints and why they exist
  - Identify potential dependency conflicts or outdated packages
  - Check for security vulnerabilities in dependencies

- [ ] **1.1.3** Create high-level architecture diagram
  - Map relationships between counting/analysis/mapping modules
  - Identify shared utilities and common code
  - Document data flow between modules
  - Identify tight coupling vs loose coupling

**Deliverables**:
- `docs/ARCHITECTURE.md` - High-level system overview
- `docs/DEPENDENCIES.md` - Dependency analysis and rationale
- `docs/DIRECTORY_STRUCTURE.md` - Complete codebase map

---

### **1.2 Module Deep Dive: Counting**
**Duration**: 2-3 days
**Dependencies**: 1.1 complete
**Priority**: HIGH

#### Tasks:
- [ ] **1.2.1** Read and document `src/counting/__main__.py`
  - CLI interface design
  - Command structure (Typer commands)
  - Parameter validation approach

- [ ] **1.2.2** Read and document `src/counting/run_counting.py`
  - Orchestration logic
  - Error handling patterns
  - Temporary file management
  - Output file conventions

- [ ] **1.2.3** Read and document `src/counting/count_alleles.py`
  - Core counting algorithm
  - BAM file reading strategy (pysam usage)
  - VCF parsing approach
  - Memory management patterns
  - Performance characteristics

- [ ] **1.2.4** Read and document `src/counting/filter_variant_data.py`
  - Filtering logic and criteria
  - Polars/Pandas usage patterns
  - Edge cases handled

- [ ] **1.2.5** Read and document `src/counting/parse_gene_data.py`
  - GTF/GFF3 parsing approach
  - Gene annotation data structures
  - Coordinate system handling (0-based vs 1-based)

- [ ] **1.2.6** Read and document single-cell variants
  - `run_counting_sc.py` - Single-cell orchestration
  - `count_alleles_sc.py` - Cell barcode handling
  - H5AD output format
  - AnnData structure

- [ ] **1.2.7** Identify counting module issues
  - Code smells (long functions, duplicate code)
  - Missing error handling
  - Performance bottlenecks
  - Security concerns (file path injection, etc.)
  - Missing input validation

**Deliverables**:
- `docs/modules/COUNTING_MODULE.md` - Complete module documentation
- `docs/modules/COUNTING_ISSUES.md` - Technical debt inventory
- Data flow diagrams for bulk and single-cell counting

---

### **1.3 Module Deep Dive: Analysis**
**Duration**: 2-3 days
**Dependencies**: 1.1 complete (can run parallel with 1.2)
**Priority**: HIGH

#### Tasks:
- [ ] **1.3.1** Read and document `src/analysis/__main__.py`
  - CLI commands structure
  - Parameter handling

- [ ] **1.3.2** Read and document `src/analysis/run_analysis.py`
  - Pipeline orchestration
  - Input validation
  - Output generation

- [ ] **1.3.3** Read and document `src/analysis/as_analysis.py` **[CRITICAL]**
  - Beta-binomial model implementation
  - Statistical optimization approach (scipy.optimize)
  - Dispersion estimation logic
  - Phased vs unphased genotype handling
  - Mathematical correctness verification

- [ ] **1.3.4** Read and document single-cell analysis
  - `run_analysis_sc.py` - Per-celltype analysis
  - `as_analysis_sc.py` - Statistical methods for scRNA
  - `compare_ai.py` - Differential allelic imbalance

- [ ] **1.3.5** Read and document support files
  - `filter_data.py` - Data filtering logic
  - `count_alleles.py` - Analysis-specific counting

- [ ] **1.3.6** Verify statistical correctness
  - Review beta-binomial formulation
  - Check optimization bounds and constraints
  - Validate p-value calculations
  - Review Z-score filtering logic

- [ ] **1.3.7** Identify analysis module issues
  - Code quality problems
  - Statistical bugs or edge cases
  - Performance issues with large datasets
  - Missing validations

**Deliverables**:
- `docs/modules/ANALYSIS_MODULE.md` - Complete module documentation
- `docs/modules/STATISTICAL_METHODS.md` - Mathematical methods documentation
- `docs/modules/ANALYSIS_ISSUES.md` - Technical debt inventory

---

### **1.4 Module Deep Dive: Mapping**
**Duration**: 2-3 days
**Dependencies**: 1.1 complete (can run parallel with 1.2, 1.3)
**Priority**: HIGH

#### Tasks:
- [ ] **1.4.1** Read and document `src/mapping/__main__.py`
  - CLI structure
  - Command flow

- [ ] **1.4.2** Read and document `src/mapping/run_mapping.py`
  - 3-step mapping pipeline orchestration
  - Temporary file handling
  - Integration with external aligners

- [ ] **1.4.3** Read and document `src/mapping/make_remap_reads.py`
  - Allele swapping algorithm
  - BAM read modification
  - SNP intersection logic

- [ ] **1.4.4** Read and document `src/mapping/filter_remap_reads.py`
  - Remapping validation logic
  - Position comparison approach
  - Read filtering criteria

- [ ] **1.4.5** Read and document data management
  - `wasp_data_files.py` - File I/O and metadata
  - `intersect_variant_data.py` - VCF/BAM intersection
  - `remap_utils.py` - Utility functions

- [ ] **1.4.6** Identify mapping module issues
  - Code quality problems
  - Edge cases in allele swapping
  - Performance with large BAMs
  - Error handling gaps

**Deliverables**:
- `docs/modules/MAPPING_MODULE.md` - Complete module documentation
- `docs/modules/MAPPING_ISSUES.md` - Technical debt inventory
- Sequence diagrams for 3-step mapping process

---

### **1.5 Cross-Cutting Concerns Analysis**
**Duration**: 1 day
**Dependencies**: 1.2, 1.3, 1.4 complete
**Priority**: MEDIUM

#### Tasks:
- [ ] **1.5.1** Identify code duplication across modules
  - Common patterns that should be utilities
  - Duplicated logic that should be shared
  - Inconsistent implementations of same concept

- [ ] **1.5.2** Analyze error handling patterns
  - Consistency of error handling
  - Quality of error messages
  - Exception types used
  - Recovery strategies

- [ ] **1.5.3** Review logging approach
  - Logging consistency
  - Log levels used appropriately
  - Sensitive data in logs

- [ ] **1.5.4** Analyze configuration management
  - Hard-coded values that should be configurable
  - Configuration file approach
  - Environment variable usage

- [ ] **1.5.5** Review data validation
  - Input validation completeness
  - File format validation
  - Parameter range checking
  - Defensive programming practices

**Deliverables**:
- `docs/CROSS_CUTTING_ISSUES.md` - Common problems and patterns
- `docs/REFACTORING_OPPORTUNITIES.md` - Consolidation opportunities

---

### **1.6 Test Coverage Analysis**
**Duration**: 1 day
**Dependencies**: 1.2, 1.3, 1.4 complete
**Priority**: HIGH

#### Tasks:
- [ ] **1.6.1** Inventory existing tests (if any)
  - What's currently tested
  - Test quality assessment
  - Test data availability

- [ ] **1.6.2** Identify critical paths requiring tests
  - Core algorithms (counting, statistical analysis, mapping)
  - Data parsing (VCF, BAM, GTF)
  - File I/O operations
  - Error conditions

- [ ] **1.6.3** Plan test strategy
  - Unit test scope and approach
  - Integration test scenarios
  - End-to-end test workflows
  - Test data requirements
  - Mock/fixture strategy

- [ ] **1.6.4** Identify test data needs
  - Can we use existing test bundle?
  - Need for synthetic test data
  - Edge case test data (empty files, malformed data)

**Deliverables**:
- `docs/testing/TEST_STRATEGY.md` - Comprehensive test plan
- `docs/testing/TEST_DATA_REQUIREMENTS.md` - Test data specifications
- `docs/testing/CURRENT_COVERAGE.md` - Existing test inventory

---

### **1.7 Security & Performance Audit**
**Duration**: 1 day
**Dependencies**: 1.2, 1.3, 1.4 complete
**Priority**: CRITICAL

#### Tasks:
- [ ] **1.7.1** Security audit
  - Command injection vulnerabilities (Bash execution, external tools)
  - Path traversal vulnerabilities
  - File permission issues
  - Temporary file security
  - Input sanitization gaps

- [ ] **1.7.2** Performance analysis
  - Memory usage patterns (large file handling)
  - CPU-intensive operations
  - I/O bottlenecks
  - Opportunities for parallelization
  - Opportunities for lazy loading/streaming

- [ ] **1.7.3** Resource management review
  - File handle leaks
  - Memory leaks
  - Proper cleanup in error paths
  - Context manager usage

**Deliverables**:
- `docs/security/SECURITY_AUDIT.md` - Security findings and recommendations
- `docs/performance/PERFORMANCE_ANALYSIS.md` - Performance bottlenecks and optimizations

---

### **1.8 Phase 1 Synthesis & Prioritization**
**Duration**: 1 day
**Dependencies**: ALL Phase 1 tasks complete
**Priority**: CRITICAL

#### Tasks:
- [ ] **1.8.1** Consolidate all findings
  - Merge all issue documents
  - Remove duplicates
  - Categorize by severity and impact

- [ ] **1.8.2** Prioritize issues for Phase 2
  - Critical: Security, correctness bugs
  - High: Performance issues, major code smells
  - Medium: Code quality, documentation
  - Low: Nice-to-haves, minor refactoring

- [ ] **1.8.3** Create Phase 2 work breakdown
  - Estimate effort for each issue
  - Identify dependencies
  - Create sprint/iteration plan

- [ ] **1.8.4** Get stakeholder alignment
  - Review findings with team/users
  - Validate priorities
  - Adjust plan based on feedback

**Deliverables**:
- `docs/PHASE1_FINDINGS.md` - Executive summary of all findings
- `docs/PHASE2_ROADMAP.md` - Prioritized work plan for Phase 2

---

## 🔧 PHASE 2: CLEAN & REFACTOR

### **2.1 Critical Security & Correctness Fixes**
**Duration**: 2-3 days
**Dependencies**: 1.8 complete
**Priority**: CRITICAL - Must be done first

#### Tasks:
- [ ] **2.1.1** Fix security vulnerabilities
  - Address command injection risks
  - Fix path traversal issues
  - Secure temporary file handling
  - Add input sanitization

- [ ] **2.1.2** Fix correctness bugs
  - Statistical calculation errors
  - Off-by-one errors in coordinate systems
  - Data corruption bugs
  - Logic errors in algorithms

- [ ] **2.1.3** Add critical input validation
  - File format validation
  - Parameter range validation
  - Dependency validation (required files exist)

- [ ] **2.1.4** Write tests for fixed bugs
  - Regression tests for each bug
  - Edge case coverage

**Deliverables**:
- All critical bugs fixed and tested
- `CHANGELOG.md` updated with fixes
- Security audit sign-off

---

### **2.2 Test Suite Foundation**
**Duration**: 3-5 days
**Dependencies**: 2.1 complete
**Priority**: CRITICAL - Needed before refactoring

#### Tasks:
- [ ] **2.2.1** Set up testing infrastructure
  - Configure pytest
  - Set up test directory structure
  - Create conftest.py with common fixtures
  - Set up CI/CD for automated testing

- [ ] **2.2.2** Create test data fixtures
  - Extract/create minimal test BAM files
  - Create minimal test VCF files
  - Create test GTF/GFF3 files
  - Create synthetic edge case data

- [ ] **2.2.3** Write unit tests for utilities
  - File parsing functions
  - Data validation functions
  - Coordinate conversion functions
  - Helper utilities

- [ ] **2.2.4** Write unit tests for core algorithms
  - Counting logic tests
  - Statistical analysis tests (with known outputs)
  - Allele swapping logic tests
  - Filtering logic tests

- [ ] **2.2.5** Write integration tests
  - End-to-end counting workflow
  - End-to-end analysis workflow
  - End-to-end mapping workflow
  - Cross-module integration points

- [ ] **2.2.6** Set up coverage reporting
  - Configure coverage.py
  - Set coverage targets (aim for >80%)
  - Integrate with CI/CD

**Deliverables**:
- Comprehensive test suite
- CI/CD pipeline running tests
- Coverage report showing >80% coverage
- `docs/testing/RUNNING_TESTS.md` - Test execution guide

---

### **2.3 Code Quality Improvements - Counting Module**
**Duration**: 2-3 days
**Dependencies**: 2.2 complete (tests provide safety net)
**Priority**: HIGH

#### Tasks:
- [ ] **2.3.1** Extract duplicated code into utilities
  - Common file I/O patterns
  - Repeated validation logic
  - Shared data transformations

- [ ] **2.3.2** Simplify complex functions
  - Break down long functions (>50 lines)
  - Extract nested logic
  - Reduce cyclomatic complexity

- [ ] **2.3.3** Improve naming conventions
  - Rename unclear variables
  - Use descriptive function names
  - Follow PEP 8 conventions

- [ ] **2.3.4** Add type hints
  - Function signatures
  - Class attributes
  - Return types

- [ ] **2.3.5** Add/improve docstrings
  - Module-level docstrings
  - Function/class docstrings
  - Follow Google/NumPy style
  - Include parameter descriptions and examples

- [ ] **2.3.6** Remove dead code
  - Unused imports
  - Commented-out code
  - Unused functions/variables

**Deliverables**:
- Refactored counting module
- All tests still passing
- Type checking passes (mypy)
- Code quality metrics improved (pylint/flake8)

---

### **2.4 Code Quality Improvements - Analysis Module**
**Duration**: 2-3 days
**Dependencies**: 2.2 complete, 2.3 complete (can learn from counting refactor)
**Priority**: HIGH

#### Tasks:
- [ ] **2.4.1** Refactor statistical analysis code
  - Extract beta-binomial model into clear functions
  - Simplify optimization logic
  - Add numerical stability checks

- [ ] **2.4.2** Apply same improvements as 2.3.1-2.3.6
  - Extract duplicated code
  - Simplify complex functions
  - Improve naming
  - Add type hints
  - Add/improve docstrings
  - Remove dead code

- [ ] **2.4.3** Improve error handling in analysis
  - Handle optimization failures gracefully
  - Validate statistical assumptions
  - Provide clear error messages for invalid inputs

**Deliverables**:
- Refactored analysis module
- All tests still passing
- Statistical methods clearly documented

---

### **2.5 Code Quality Improvements - Mapping Module**
**Duration**: 2-3 days
**Dependencies**: 2.2 complete, 2.3, 2.4 complete
**Priority**: HIGH

#### Tasks:
- [ ] **2.5.1** Refactor allele swapping logic
  - Clarify algorithm steps
  - Improve edge case handling
  - Add validation

- [ ] **2.5.2** Apply same improvements as 2.3.1-2.3.6
  - Extract duplicated code
  - Simplify complex functions
  - Improve naming
  - Add type hints
  - Add/improve docstrings
  - Remove dead code

**Deliverables**:
- Refactored mapping module
- All tests still passing

---

### **2.6 Error Handling & Logging Standardization**
**Duration**: 2 days
**Dependencies**: 2.3, 2.4, 2.5 complete
**Priority**: MEDIUM

#### Tasks:
- [ ] **2.6.1** Standardize error handling
  - Define custom exception hierarchy
  - Consistent try/except patterns
  - Proper error propagation

- [ ] **2.6.2** Improve error messages
  - Make messages actionable
  - Include context (file names, line numbers)
  - Suggest fixes where possible

- [ ] **2.6.3** Standardize logging
  - Consistent log format
  - Appropriate log levels
  - Structured logging (JSON option)
  - Remove sensitive data from logs

- [ ] **2.6.4** Add progress indicators
  - Progress bars for long operations
  - Status updates for multi-step workflows
  - Estimated time remaining where possible

**Deliverables**:
- Standardized error handling across codebase
- Consistent logging approach
- Better user experience during execution

---

### **2.7 Configuration & Flexibility Improvements**
**Duration**: 1-2 days
**Dependencies**: 2.3, 2.4, 2.5 complete
**Priority**: MEDIUM

#### Tasks:
- [ ] **2.7.1** Extract hard-coded values
  - Create configuration file(s)
  - Use environment variables appropriately
  - Document all configuration options

- [ ] **2.7.2** Add configuration validation
  - Validate config on startup
  - Provide helpful errors for invalid configs

- [ ] **2.7.3** Create config examples
  - Default configuration
  - Example configurations for common use cases

**Deliverables**:
- Configuration system implemented
- `docs/CONFIGURATION.md` - Configuration guide
- Example config files

---

### **2.8 Documentation Finalization**
**Duration**: 2-3 days
**Dependencies**: All Phase 2 refactoring complete
**Priority**: HIGH

#### Tasks:
- [ ] **2.8.1** Update README.md
  - Installation instructions
  - Quick start guide
  - Link to detailed docs

- [ ] **2.8.2** Create user documentation
  - Tutorial for each module
  - Example workflows
  - Common use cases
  - Troubleshooting guide

- [ ] **2.8.3** Create developer documentation
  - Contributing guide
  - Code style guide
  - Architecture overview
  - Module interaction diagrams

- [ ] **2.8.4** API documentation
  - Generate API docs from docstrings (Sphinx)
  - Host docs (Read the Docs or GitHub Pages)

- [ ] **2.8.5** Create CHANGELOG
  - Document all changes from Phase 1 & 2
  - Version appropriately (semantic versioning)

**Deliverables**:
- Complete user documentation
- Complete developer documentation
- Generated API documentation
- Updated README and CHANGELOG

---

### **2.9 Performance Optimization**
**Duration**: 2-3 days
**Dependencies**: 2.3, 2.4, 2.5 complete (correctness first)
**Priority**: MEDIUM (can be done in parallel with 2.6-2.8)

#### Tasks:
- [ ] **2.9.1** Profile critical paths
  - Use cProfile on representative workloads
  - Identify bottlenecks

- [ ] **2.9.2** Optimize I/O operations
  - Implement streaming where possible
  - Use efficient file readers
  - Reduce disk I/O

- [ ] **2.9.3** Optimize data structures
  - Use appropriate data structures (Polars vs Pandas)
  - Reduce memory footprint
  - Implement chunking for large files

- [ ] **2.9.4** Add parallelization where appropriate
  - Identify parallelizable operations
  - Use multiprocessing/threading carefully
  - Benchmark improvements

- [ ] **2.9.5** Optimize algorithms
  - Improve algorithmic complexity where possible
  - Use vectorized operations
  - Leverage NumPy/Polars optimizations

**Deliverables**:
- Performance benchmarks (before/after)
- Optimized code with tests still passing
- `docs/performance/BENCHMARKS.md` - Performance metrics

---

### **2.10 Phase 2 Validation & Release Prep**
**Duration**: 1-2 days
**Dependencies**: ALL Phase 2 tasks complete
**Priority**: CRITICAL

#### Tasks:
- [ ] **2.10.1** Comprehensive testing
  - Run full test suite
  - Manual testing of key workflows
  - Test on different platforms if applicable

- [ ] **2.10.2** Code review
  - Internal code review
  - External review if possible
  - Address all findings

- [ ] **2.10.3** Documentation review
  - Ensure all docs are accurate
  - Check for broken links
  - Verify examples work

- [ ] **2.10.4** Prepare release
  - Tag version in git
  - Update version numbers
  - Create release notes

- [ ] **2.10.5** Acceptance testing
  - Test with real user workflows
  - Verify all use cases still work
  - Get user feedback

**Deliverables**:
- Release-ready codebase
- All tests passing
- All documentation complete
- Ready for Phase 3 (feature additions)

---

## 📊 Success Metrics

### Phase 1 Complete When:
- [ ] All modules documented
- [ ] All technical debt identified and prioritized
- [ ] Test strategy defined
- [ ] Security audit complete
- [ ] Phase 2 roadmap approved

### Phase 2 Complete When:
- [ ] 0 critical security issues
- [ ] >80% test coverage
- [ ] All code has type hints
- [ ] All code has docstrings
- [ ] Code quality metrics improved (pylint score >8.0)
- [ ] Performance benchmarks meet targets
- [ ] All documentation complete
- [ ] User acceptance testing passed

---

## 🚀 Estimated Timeline

| Phase | Tasks | Duration | Dependencies |
|-------|-------|----------|--------------|
| **Phase 1** | 1.1 - 1.8 | **10-14 days** | None |
| **Phase 2** | 2.1 - 2.10 | **20-28 days** | Phase 1 complete |
| **TOTAL** | | **30-42 days** (~6-8 weeks) | |

---

## ⚠️ Risk Mitigation

### Risks & Mitigation Strategies:

1. **Breaking existing functionality during refactor**
   - Mitigation: Write comprehensive tests BEFORE refactoring (2.2)
   - Run tests after every change

2. **Scope creep during exploration**
   - Mitigation: Stick to documentation in Phase 1, no fixes
   - Save all improvements for Phase 2

3. **Statistical algorithm changes introducing bugs**
   - Mitigation: Validate against known test cases
   - Get statistical review from domain expert
   - Maintain backward compatibility

4. **Documentation becoming outdated**
   - Mitigation: Update docs as part of each task
   - Include doc updates in PR reviews

5. **Performance regressions from refactoring**
   - Mitigation: Establish baseline benchmarks early
   - Monitor performance in CI/CD
   - Profile before and after major changes

---

## 🔄 Continuous Practices

Throughout Phases 1 & 2:

- [ ] Use feature branches for all changes
- [ ] Require code review before merging
- [ ] Run automated tests on every commit
- [ ] Update documentation with code changes
- [ ] Keep stakeholders informed of progress
- [ ] Hold regular retrospectives to improve process

---

## 📝 Notes

- This plan assumes 1-2 developers working full-time
- Durations are estimates; adjust based on actual findings
- Some tasks can be parallelized to reduce total time
- Priority markers guide what to do first if resources are limited
- After Phase 2, codebase will be ready for PGEN and other feature additions

---

**Next Step**: Get approval for this plan, then begin with Task 1.1.1
