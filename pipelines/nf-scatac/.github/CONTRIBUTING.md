# Contributing to nf-scatac

## Getting help

For questions, bugs, or feature requests, please open an issue on [GitHub](https://github.com/mcvickerlab/WASP2/issues).

## Development workflow

1. Fork the repository
2. Create a feature branch from `dev`
3. Make your changes
4. Run `nf-core pipelines lint` to verify compliance
5. Submit a pull request to `dev`

## Code style

- Follow nf-core module conventions for new modules
- Use `tuple val(meta), path(...)` for all process inputs/outputs
- Include `stub:` blocks in all processes
- Add `versions.yml` output to all processes
- Write `meta.yml` documentation for new modules

## Testing

Run the test profile before submitting changes:

```bash
nextflow run main.nf -profile test,docker --outdir test_results
```
