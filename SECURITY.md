# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.3.x   | :white_check_mark: |
| 1.2.x   | :white_check_mark: |
| < 1.2   | :x:                |

## Reporting a Vulnerability

We take security vulnerabilities seriously. If you discover a security issue, please report it responsibly.

### How to Report

1. **Do NOT create a public GitHub issue** for security vulnerabilities
2. Use [GitHub's private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing/privately-reporting-a-security-vulnerability) or email the maintainers directly
3. Include:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Any suggested fixes (optional)

### What to Expect

- **Acknowledgment**: Within 48 hours of your report
- **Initial Assessment**: Within 1 week
- **Resolution Timeline**: Depends on severity
  - Critical: 1-2 weeks
  - High: 2-4 weeks
  - Medium/Low: Next release cycle

### Disclosure Policy

- We follow coordinated disclosure practices
- We will credit reporters in release notes (unless anonymity is requested)
- Please allow us reasonable time to address issues before public disclosure

## Security Measures

This project implements multiple security scanning tools:

### Dependency Scanning
- **Dependabot**: Automatic security updates for Python (pip), Rust (cargo), and GitHub Actions
- **pip-audit**: Python dependency vulnerability scanning in CI
- **cargo-audit**: Rust dependency vulnerability scanning in CI

### Static Analysis
- **Bandit**: Python security linter (configured in pre-commit and CI)
- **CodeQL**: GitHub's advanced static analysis for Python
- **Ruff**: Fast Python linter and formatter

### Secret Detection
- **Gitleaks**: Pre-commit hook for detecting secrets and credentials
- **detect-private-key**: Pre-commit hook for private key detection

### Container Security
- Multi-stage Docker builds with non-root user
- Minimal base images (python:3.11-slim)

## Security Best Practices for Contributors

1. **Never commit secrets** - Use environment variables or secret management
2. **Keep dependencies updated** - Review and merge Dependabot PRs promptly
3. **Run pre-commit hooks** - Ensures security checks pass locally
4. **Review security alerts** - Check GitHub Security tab regularly

## Security Contacts

For security-related inquiries, contact the project maintainers through GitHub's private vulnerability reporting feature or the repository's security advisories.
