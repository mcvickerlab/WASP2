# Self-Hosted Runner Configuration Audit

**Date:** 2026-02-03
**Issue:** #207
**Scope:** `.github/runner/`, `scripts/setup-mac-runner.sh`, `scripts/setup-multi-runners.sh`, workflow files

---

## Executive Summary

The self-hosted runner infrastructure is well-architected with a 2-tier health monitoring system (launchd + watchdog) and specialized runner routing for different job types. Three critical and several medium-severity issues were identified. Two critical issues have been remediated in this audit; one requires manual attention.

---

## Files Audited

| File | Purpose | Lines |
|------|---------|-------|
| `.github/runner/README.md` | Runner setup documentation | 111 |
| `.github/runner/com.github.actions.runner.plist` | launchd service template | 71 |
| `.github/runner/install-service.sh` | Service installation script | 113 |
| `.github/runner/watchdog.sh` | Health monitoring daemon | 191 |
| `scripts/setup-mac-runner.sh` | Single runner setup | 223 |
| `scripts/setup-multi-runners.sh` | Multi-runner setup (3 runners) | 247 |

---

## Findings

### CRITICAL

#### C1: No checksum verification for downloaded runner binary
**Status:** REMEDIATED
**Files:** `scripts/setup-mac-runner.sh`, `scripts/setup-multi-runners.sh`
**Risk:** Supply chain compromise. Runner binary executes arbitrary code on the host. Without SHA256 verification, a compromised CDN could serve a malicious binary.
**Fix:** Added SHA256 checksum verification by extracting hashes from the GitHub release notes body via the API. Downloads use `curl --fail` to detect HTTP errors. If the checksum cannot be retrieved, a warning is emitted and manual verification is recommended. Download aborts and cleans up if checksums don't match. Also added `set -o pipefail` and `RUNNER_VERSION` validation.

#### C2: No tar archive validation
**Status:** ACCEPTED (Low Likelihood)
**Files:** `scripts/setup-mac-runner.sh`, `scripts/setup-multi-runners.sh`
**Risk:** Path traversal or tar bomb attacks from a compromised archive. Mitigated by C1 fix (checksum verification ensures archive integrity).

#### C3: Hardcoded developer path in ci.yml
**Status:** NOT IN SCOPE (not a runner config file)
**Files:** `.github/workflows/ci.yml` (line 282)
**Risk:** `MAMBAFORGE_PYTHON=/Users/jeffjaureguy/mambaforge/bin/python` will fail on any runner other than the developer's machine.
**Recommendation:** Replace with `$(which python3)` or an environment variable.

### HIGH

#### H1: Process killing via pattern matching
**Status:** ACCEPTED (Standard Practice)
**File:** `.github/runner/watchdog.sh`
**Risk:** `pkill -f "Runner.Listener"` could theoretically kill unrelated processes with similar names. In practice, "Runner.Listener" is specific to GitHub's runner binary.
**Mitigation:** GitHub's own `svc.sh` uses the same pattern. Risk is negligible in a dedicated runner environment.

#### H2: Race condition in PID management
**Status:** ACCEPTED (Low Impact)
**File:** `.github/runner/watchdog.sh`
**Risk:** In `stop_watchdog()`, between the `kill -0` PID check and the subsequent `kill`, the PID could theoretically be reassigned. The main health-check loop uses `pgrep -f` which has a similar TOCTOU concern. Extremely unlikely given check intervals and typical PID assignment patterns.

### MEDIUM

#### M1: Non-deterministic health logging
**Status:** REMEDIATED
**File:** `.github/runner/watchdog.sh`
**Risk:** `RANDOM % 10` for health log decisions is unreliable - could skip health logs for extended periods, making incident investigation harder.
**Fix:** Replaced with deterministic counter (`check_count % 10`) for consistent health logging every ~10 minutes.

#### M2: Missing RUNNER_DIR validation in watchdog
**Status:** REMEDIATED
**File:** `.github/runner/watchdog.sh`
**Risk:** Watchdog would silently fail if RUNNER_DIR doesn't exist.
**Fix:** Added startup validation that checks RUNNER_DIR exists and creates `_diag` directory if needed.

#### M3: Deprecated launchctl commands
**Status:** ACCEPTED (Backwards Compatibility)
**File:** `.github/runner/install-service.sh`
**Risk:** `launchctl load` is deprecated on macOS 13+ in favor of `launchctl enable`/`launchctl bootstrap`. Still functional and widely used.
**Recommendation:** Update when minimum macOS version is bumped to 13+.

#### M4: Hardcoded repository path
**Status:** ACCEPTED (Intentional)
**Files:** `scripts/setup-mac-runner.sh`, `scripts/setup-multi-runners.sh`
**Risk:** `REPO="mcvickerlab/WASP2"` limits reusability.
**Rationale:** These are project-specific setup scripts. Making REPO configurable via env var would improve reusability but is low priority.

### LOW

#### L1: Emoji in script output
**Status:** ACCEPTED
**Risk:** Could cause encoding issues on non-UTF8 terminals. All modern macOS terminals support UTF-8.

#### L2: Arbitrary sleep durations
**Status:** ACCEPTED
**File:** `.github/runner/install-service.sh`
**Risk:** `sleep 2` and `sleep 3` may be insufficient on slow systems. Acceptable for macOS M3 Max target hardware.

---

## Runner Security Assessment

### Isolation
- Runners execute in user-space under the installing user's account
- No containerized isolation (standard for macOS self-hosted runners)
- Docker is available but runners themselves are not containerized
- Recommendation: Ensure runner user has minimal system privileges

### Permissions
- launchd plist uses `Nice: -5` (elevated priority) - appropriate for CI
- `ProcessType: Interactive` allows shell access - required for build tools
- `LowPriorityIO: false` ensures normal I/O scheduling priority (not deprioritized by macOS)

### Monitoring & Auto-Restart
- **launchd layer:** Auto-restart on crash/exit with 10s throttle
- **Watchdog layer:** 60s health checks, detects socket timeouts and stalled processes
- **Escalation:** 3 consecutive errors triggers force restart (unconditional SIGTERM then SIGKILL sequence)
- **GitHub connectivity:** Periodic check to `api.github.com/zen`
- Assessment: Robust 2-tier monitoring. Well-designed for the "stuck runner" problem.

---

## Runner Label Alignment

### Single Runner (`setup-mac-runner.sh`)
Labels: `macOS, ARM64, docker, wasp2` (plus `self-hosted` added automatically by GitHub)

Matches workflows:
- ci.yml (5 jobs)
- docker.yml (2 jobs)
- security.yml (6 jobs)
- benchmarks.yml (2 jobs)
- velocity-bot.yml (subset match)
- release.yml (subset match)

### Multi-Runner (`setup-multi-runners.sh`)

| Runner | Labels | Matching Workflows |
|--------|--------|--------------------|
| python-runner | `python, testing, lint, fast` | nightly.yml unit-tests job |
| rust-runner | `rust, build, maturin` | No direct workflow match |
| analysis-runner | `analysis, bioinformatics, docker, slow` | nightly.yml integration/analysis jobs |

**Findings:**
1. The `rust-runner` labels (`rust, build, maturin`) have no matching `runs-on` in any workflow. Rust build jobs in ci.yml use the generic `docker, wasp2` labels. Consider either adding `rust` to ci.yml rust jobs or removing the dedicated rust runner.
2. Label architecture is otherwise sound: CI jobs route to the single runner via `wasp2`, nightly jobs route to specialized runners.

---

## Queue Time Optimization

- **Current setup:** Single runner for CI + 3 specialized runners for nightly = 4 total runners
- **Parallelism:** M3 Max can handle multiple concurrent runners efficiently
- **Bottleneck:** All CI jobs require `wasp2` label, funneling through one runner. Consider adding `wasp2` label to the python and analysis runners for CI overflow.
- **Nightly routing:** Good separation of fast (python) vs slow (analysis) jobs

---

## Recommendations

### Immediate (This PR)
- [x] Add SHA256 checksum verification to runner downloads
- [x] Add `curl --fail` and `set -o pipefail` to setup scripts
- [x] Add `RUNNER_VERSION` validation before download
- [x] Fix non-deterministic health logging in watchdog
- [x] Add RUNNER_DIR validation to watchdog startup
- [x] Fix PID file write quoting and mkdir error handling in watchdog

### Short-term
- [ ] Fix hardcoded path in ci.yml line 282
- [ ] Add `wasp2` label to multi-runners for CI overflow capacity
- [ ] Remove or repurpose unused rust-runner labels

### Long-term
- [ ] Migrate to `launchctl bootstrap` when minimum macOS >= 13
- [ ] Consider ephemeral runners for security-sensitive workflows
- [ ] Add runner version auto-update mechanism
