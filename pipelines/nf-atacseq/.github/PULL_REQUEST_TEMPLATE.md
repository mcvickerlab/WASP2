## PR checklist

- [ ] This comment contains a description of changes with context.
- [ ] Tests pass locally with `nf-test test` using the `test,docker` profile.
- [ ] If you have added new modules/processes, they include `nf-test` tests.
- [ ] If applicable, new parameters are documented in `nextflow_schema.json`.
- [ ] Pipeline runs successfully with `test` and `test_stub` profiles.
- [ ] `CHANGELOG.md` is updated with noteworthy changes.

## Description

<!-- Please include a summary of the change and which issue is fixed. -->

Fixes # (issue)

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## How has this been tested?

<!-- Please describe the tests you ran. Provide instructions so we can reproduce. -->

- [ ] `nf-test test --profile test,docker`
- [ ] `nf-test test --profile test_stub,docker`
- [ ] Manual pipeline run with real data
