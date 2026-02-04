# AGENTS.md — Velocity Bot Configuration

This document defines the behavior, permissions, and guardrails for the automated velocity-bot system.

## Overview

The velocity bot automates issue implementation using Claude Code. It triggers on:
- **Label trigger:** Adding the `velocity` label to an issue
- **Command trigger:** Commenting `/implement` on an issue

Only repository OWNER, MEMBER, or COLLABORATOR can trigger the bot.

## Architecture

```
Issue trigger → check-trigger (ubuntu-latest) → implement (self-hosted, macOS, ARM64)
```

- **check-trigger**: Lightweight event filter, auth check, input sanitization
- **implement**: Claude Code execution, commit, PR creation

## Permissions

### Allowed Tools
- File operations: `Read`, `Edit`, `Write`, `Grep`, `Glob`
- Git (read-only): `git diff`, `git log`, `git status`, `git show`
- Build/test: `pytest`, `python`, `maturin`, `cargo`, `ruff`
- Filesystem: `ls`

### Denied Tools
- Network: `curl`, `wget`, `ssh`, `scp`, `WebFetch`
- Destructive: `rm -rf`, `rm -r`, `chmod 777`, `sudo`
- Secrets: `.env*`, `~/.ssh/*`, `~/.aws/*`, `~/.gnupg/*`

## Guardrails

| Guardrail | Value |
|-----------|-------|
| Execution timeout | 30 minutes |
| Max Claude turns | 50 |
| Model | Claude Sonnet |
| Auth gating | OWNER / MEMBER / COLLABORATOR |
| Anti-loop | `github.actor != 'claude[bot]'` |
| Input sanitization | Strip HTML comments, invisible Unicode, 4000 char limit |
| Permission model | Explicit allowlist (not --dangerously-skip-permissions) |

## Labels

| Label | Purpose |
|-------|---------|
| `velocity` | Trigger automation |
| `bot:in-progress` | Bot is actively working |
| `bot:pr-ready` | PR created, awaiting review |
| `bot:failed` | Execution failed |
| `bot:needs-help` | Needs human input |

## Scope Limits

- Bot makes conservative, focused changes
- Does not refactor unrelated code
- Does not push to main directly — always creates a PR
- All changes require human review before merge

## Escalation

1. Bot posts failure details to the issue
2. `bot:needs-help` label signals human intervention needed
3. Check Actions run logs for debugging

## Break Glass

To disable the bot:
1. Remove the `velocity` label from the issue
2. Delete or disable `.github/workflows/velocity-bot.yml`
3. Revoke the `ANTHROPIC_API_KEY` secret

## Monitoring

Weekly dashboard updates posted to a pinned "Velocity Bot Dashboard" issue.
Metrics tracked: issues processed, PRs created, merge rate, failure categories.
