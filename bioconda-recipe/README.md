# WASP2 Bioconda Recipe

This directory contains the Bioconda recipe for WASP2.

## Submission Process

After the package is published to PyPI:

1. **Fork bioconda-recipes**
   ```bash
   gh repo fork bioconda/bioconda-recipes --clone
   cd bioconda-recipes
   ```

2. **Create recipe directory**
   ```bash
   mkdir -p recipes/wasp2
   cp /path/to/wasp2/bioconda-recipe/meta.yaml recipes/wasp2/
   cp /path/to/wasp2/bioconda-recipe/build.sh recipes/wasp2/
   ```

3. **Update sha256 hash**
   Get the hash from PyPI:
   ```bash
   curl -sL https://pypi.io/packages/source/w/wasp2/wasp2-1.3.0.tar.gz | sha256sum
   ```
   Update the `sha256:` field in `meta.yaml`.

4. **Test locally with bioconda-utils**
   ```bash
   # Install bioconda-utils
   conda create -n bioconda -c conda-forge -c bioconda bioconda-utils
   conda activate bioconda

   # Lint the recipe
   bioconda-utils lint --packages wasp2

   # Build locally
   bioconda-utils build --packages wasp2
   ```

5. **Submit PR**
   ```bash
   git checkout -b add-wasp2
   git add recipes/wasp2/
   git commit -m "Add wasp2 recipe"
   git push origin add-wasp2
   gh pr create --title "Add wasp2 1.3.0" --body "New recipe for WASP2: allele-specific analysis of NGS data"
   ```

## Recipe Notes

- Uses `{{ compiler('rust') }}` for Rust toolchain
- Uses `cargo-bundle-licenses` to bundle Rust dependency licenses (Bioconda requirement)
- Requires htslib for rust-htslib compilation
- Skips Windows and Python <3.10

## Testing

```bash
conda create -n wasp2-test -c bioconda -c conda-forge wasp2
conda activate wasp2-test
wasp2-count --help
python -c "import wasp2_rust; print('OK')"
```
