# Preventing AI Agent File Pollution: Deep Research Report

**Date:** 2026-03-06
**Context:** WASP2 Nextflow pipelines — cleanup after nf-core lint compliance work
**Research method:** 7 specialized agents covering Anthropic CLI, git worktrees, Agent SDK, pre-commit hooks, community practices, CLAUDE.md patterns, and Nextflow-specific cleanup

---

## Executive Summary

AI coding agents (Claude Code, Cursor, Copilot, Devin) systematically create files that were never requested — placeholder images, temp files, debug scripts, documentation, and binary artifacts. **No single mechanism prevents this.** The industry consensus is **defense-in-depth**: layer multiple controls so no single failure lets junk through.

Claude Code has **no built-in cleanup mechanism**. The `SessionEnd` hook is the only lifecycle event for cleanup, and it must be configured manually. File checkpointing exists but only tracks Write/Edit tool operations, not Bash commands.

---

## The Problem We Experienced

During WASP2 nf-core compliance work, agents created:
- `-.bam` (272MB) — stray BAM from a failed pipe command
- 9 placeholder PNGs (69-113 bytes) — broken images in `assets/` and `docs/images/`
- 20 `CLAUDE.md` files — auto-generated activity logs in pipeline subdirectories
- Various `work/`, `.nextflow/`, `results_*/` directories from test runs

---

## What Anthropic Provides (and Doesn't)

### Available Mechanisms

| Mechanism | What It Does | Limitations |
|-----------|-------------|-------------|
| **SessionEnd hook** | Fires when session terminates | Must configure manually; async; no guarantee of completion |
| **WorktreeRemove hook** | Fires when worktree is destroyed | Only for worktree-isolated sessions |
| **File checkpointing** | Snapshots files before Write/Edit; can rewind | Does NOT track Bash-created files (`echo >`, `cp`, `mv`) |
| **Sandbox** (macOS Seatbelt / Linux bwrap) | Restricts where files can be written | Prevents unauthorized writes but doesn't clean up authorized ones |
| **PreToolUse hooks** | Can block Write operations via pattern matching | Known reliability issues (GitHub #16733, #4362) |
| **Worktree isolation** | Each agent gets separate working directory | Orphaned worktrees after crashes (GitHub #26725) |

### What's Missing
- No automatic file tracking or artifact registry
- No cleanup-on-exit by default
- No diff-from-baseline mechanism
- No "created files" manifest per session
- Sandbox mode bugs create stray 0-byte files (GitHub #18548)

---

## Defense-in-Depth Strategy (Recommended)

### Layer 1: `.gitignore` (Passive Defense)

Comprehensive `.gitignore` prevents accidental commits. Our WASP2 `.gitignore` now covers:

```gitignore
# Nextflow runtime
.nextflow/
.nextflow.log*
work/

# Test/results artifacts
test_data/
test_results/
test-output/
results_stub/
pipelines/*/results_*/

# Claude Code state
.claude/
**/CLAUDE.md
!./CLAUDE.md

# Nextflow reports
trace.txt
timeline.html
report.html
dag.svg
dag.dot

# Large binaries
*.bam
*.bam.bai
*.vcf.gz
*.vcf.gz.tbi
```

### Layer 2: CLAUDE.md Directives (Behavioral)

Add explicit file creation rules to project CLAUDE.md:

```markdown
## File Hygiene Rules
- NEVER create files unless absolutely necessary for the task
- NEVER create placeholder/stub files (empty PNGs, dummy data)
- NEVER create files in the repo root — use appropriate subdirectories
- ALWAYS prefer editing existing files over creating new ones
- ALWAYS clean up temp files created during debugging
- If you create a test script, delete it when done
```

**Effectiveness:** Inconsistent. GitHub issue #2901 documents systematic CLAUDE.md instruction violations. Use as one layer, not the only layer.

**Key insight from community:** "If a rule is important enough to say NEVER in CLAUDE.md, it should be enforced by a hook, not just text." ([DEV Community](https://dev.to/siddhantkcode/an-easy-way-to-stop-claude-code-from-forgetting-the-rules-h36))

**CLAUDE.md sizing:** Keep under 80 lines total. Over that, Claude starts ignoring parts. Use `@imports` for overflow. ([HumanLayer](https://www.humanlayer.dev/blog/writing-a-good-claude-md))

**Purpose-built plugin:** [claude-anti-flood](https://github.com/angeloticiano/claude-anti-flood) — 4 modular SKILLs (prevent-file-floods, prevent-skill-bloat, prevent-agent-bloat, prevent-command-bloat). Enforces "chat-first" approach: deliver analysis inline, not as files.

**nf-core stance on AI files:** "Aims to not clutter the template with AI helper files (AGENT.md, CLAUDE.md, .mcp.json)" ([nf-core blog](https://nf-co.re/blog/2026/statement-on-ai))

### Layer 3: Pre-commit Hooks (Automated Gate)

Install `pre-commit` with hooks that catch agent artifacts:

```yaml
# .pre-commit-config.yaml
repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v5.0.0
    hooks:
      - id: check-added-large-files
        args: ['--maxkb=500']
      - id: check-merge-conflict
      - id: check-symlinks
      - id: detect-private-key
      - id: end-of-file-fixer
      - id: trailing-whitespace

  # Block known agent artifact patterns
  - repo: local
    hooks:
      - id: forbid-agent-artifacts
        name: Block AI agent artifacts
        entry: "Forbidden file: likely an AI agent artifact"
        language: fail
        files: |
          (?x)^(
            ANALYSIS\.md|
            CHANGELOG\.md|
            debug_.*\.py|
            test_scratch.*|
            tmpclaude.*|
            .*\.debug\.(py|js|ts)|
            run_.*\.py
          )$
```

**Agent-Precommit** (purpose-built tool): Runs differentiated checks — fast (~3s) for humans, thorough (~5min) for agent commits. Auto-detects agent sessions via `CLAUDE_CODE`, `CURSOR_SESSION`, `AIDER_MODEL` env vars.

### Layer 4: SessionEnd Cleanup Hook (Automated)

Create a cleanup script that runs when Claude Code exits:

```json
// .claude/settings.json
{
  "hooks": {
    "SessionEnd": [
      {
        "type": "command",
        "command": "bash .claude/hooks/session-cleanup.sh"
      }
    ]
  }
}
```

```bash
#!/usr/bin/env bash
# .claude/hooks/session-cleanup.sh
# Runs on every Claude Code session exit

cd "${CLAUDE_CWD:-$(pwd)}"

# Remove known agent artifacts
find . -name "tmpclaude-*-cwd" -delete 2>/dev/null
find . -name "*.debug.py" -delete 2>/dev/null
find . -name "ANALYSIS.md" -delete 2>/dev/null

# Remove empty files (likely stubs)
find . -type f -empty -not -path './.git/*' -not -name '.gitkeep' -delete 2>/dev/null

# Remove placeholder images (< 200 bytes)
find . -name "*.png" -size -200c -not -path './.git/*' -delete 2>/dev/null

# Clean Nextflow artifacts from development runs
find . -name "trace.txt" -delete 2>/dev/null
find . -name "timeline.html" -delete 2>/dev/null
find . -name "report.html" -delete 2>/dev/null
find . -name "dag.svg" -delete 2>/dev/null
find . -name "dag.dot" -delete 2>/dev/null

# Report what was cleaned
echo "Session cleanup complete"
```

### Layer 5: Git Discipline (Process)

1. **Commit before every agent interaction** — enables instant `git checkout .` rollback
2. **Never `git add .`** — stage files individually
3. **`git clean -n` after sessions** — preview untracked files before deleting
4. **`git diff --stat` before commit** — verify no unexpected files

### Layer 6: Worktree Isolation (Structural)

For parallel agent work, use `--worktree` to isolate each agent:

```bash
claude --worktree feat-lint-fixes
```

Or configure subagents with `isolation: "worktree"` in agent frontmatter. Each worktree gets its own branch and working directory.

**Cleanup orphaned worktrees periodically:**
```bash
git worktree list          # Find orphans
git worktree prune         # Remove stale entries
```

---

## Nextflow-Specific Cleanup

### `nextflow clean` Command
```bash
nextflow clean -f                           # Clean latest run
nextflow clean -f -before $(nextflow log -q | tail -n 1)  # Clean all but latest
nextflow clean -n                           # Dry run (preview)
nextflow clean -f -keep-logs               # Preserve .command.log files
```

### Canonical nf-core `.gitignore` entries
```
.nextflow*
work/
results/
data/
testing/
null/
```

### Config-based cleanup
```groovy
// nextflow.config — auto-clean work dirs on success
cleanup = true
```

---

## Comparison: How Other AI Tools Handle Isolation

| Tool | Isolation Method | Cleanup |
|------|-----------------|---------|
| **Claude Code** | Git worktree | Auto if no changes; prompt if changes; orphans on crash |
| **Cursor** | Git worktree per parallel agent | Manual branch merge |
| **OpenAI Codex** | Container sandbox (bwrap/seccomp) | Container destroyed after task |
| **GitHub Copilot Agent** | Ephemeral GitHub Actions VM | VM destroyed after task |
| **Devin** | Full isolated VM ("Devbox") | VM ephemeral |

Key insight: Container/VM-based tools solve cleanup automatically (destroy the container). Worktree-based tools require explicit cleanup.

---

## Community-Reported Issues (with GitHub Issue References)

| Issue | Description | Link |
|-------|-------------|------|
| Temp file litter | 155 `tmpclaude-*-cwd` files accumulated | [#17720](https://github.com/anthropics/claude-code/issues/17720) |
| Sandbox 0-byte files | Sandbox mode creates stray empty files | [#18548](https://github.com/anthropics/claude-code/issues/18548) |
| CLAUDE.md violations | Systematic non-compliance with file rules | [#2901](https://github.com/anthropics/claude-code/issues/2901) |
| PreToolUse unreliable | Hooks "never invoked", files created anyway | [#16733](https://github.com/anthropics/claude-code/issues/16733) |
| Orphaned worktrees | Never cleaned up after crashes | [#26725](https://github.com/anthropics/claude-code/issues/26725) |
| `approve: false` ignored | Hook deny decisions not enforced | [#4362](https://github.com/anthropics/claude-code/issues/4362) |

---

## Recommended Implementation for WASP2

### Immediate (do now)
1. Add `session-cleanup.sh` hook to `.claude/hooks/`
2. Add file hygiene rules to WASP2 `CLAUDE.md`
3. Install `pre-commit` with `check-added-large-files` and `forbid-agent-artifacts`

### Short-term
4. Configure `PreToolUse` hook to block Write operations for `*.bam`, `*.png` in wrong directories
5. Add `nextflow clean` to development workflow documentation
6. Periodic `git worktree prune` in development workflow

### Long-term
7. Evaluate `agent-precommit` for differentiated human/agent commit gating
8. Consider container-based isolation for batch agent workflows

---

## Key Takeaway

> The developers who report the least friction use **defense-in-depth**: explicit rules AND hooks AND git discipline AND pre-commit checks AND scoped tasks, layered together so no single failure point lets junk through.

No single mechanism is reliable. Layer them.
