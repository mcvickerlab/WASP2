# WASP2 — Project Instructions

## What This Is
WASP2 is a mapping bias correction tool for genomic analyses. It includes:
- A Python/Rust core library (`src/wasp2/`)
- 4 Nextflow pipelines (`pipelines/nf-atacseq`, `nf-rnaseq`, `nf-scatac`, `nf-outrider`)
- Benchmarking infrastructure (`benchmarking/`)

## File Hygiene Rules
- NEVER create files unless absolutely necessary for the task
- NEVER create placeholder/stub files (empty PNGs, dummy data, skeleton configs)
- NEVER create files in the repo root — use appropriate subdirectories
- ALWAYS prefer editing existing files over creating new ones
- ALWAYS clean up temp files created during debugging before finishing
- If you create a test/debug script, delete it when done
- NEVER commit binary files (BAM, BAI, VCF, FASTQ) to git — use test fixtures in `tests/`

## Nextflow Development
- Pipelines follow nf-core conventions (modules, subworkflows, configs)
- Use `nextflow clean -f` after test runs to remove work directories
- Test profiles use small chr21 data — don't download full genomes
- Aligner: BWA (default) or Bowtie2 via `--aligner` parameter

## Code Style
- Python: ruff for linting/formatting, basedpyright for types
- Pre-commit hooks are configured — run `pre-commit install` after cloning
- Nextflow: follow nf-core module patterns (meta map, versions.yml emit)

## Git Workflow
- Feature branches → PR to `dev` → release merges `dev` → `main`
- Never commit directly to `main` or `dev`
- Stage files individually (`git add <file>`) — never `git add .`
- Run `pre-commit run --all-files` before committing
