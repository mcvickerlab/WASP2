# Bioconda Recipe for WASP2

This directory contains the Bioconda recipe template for WASP2.

## Prerequisites

Before submitting to Bioconda, ensure:

1. **PyPI package is published**: `pip install wasp2` must work
2. **Get sha256 hash**: After PyPI publish, run:
   ```bash
   pip download wasp2==1.3.0 --no-binary :all: --no-deps
   sha256sum wasp2-1.3.0.tar.gz
   ```

## Submission Steps

1. **Fork bioconda-recipes**
   ```bash
   gh repo fork bioconda/bioconda-recipes --clone
   cd bioconda-recipes
   ```

2. **Create recipe directory**
   ```bash
   mkdir -p recipes/wasp2
   cp /path/to/this/meta.yaml recipes/wasp2/
   ```

3. **Update sha256** in `meta.yaml` with actual hash from PyPI

4. **Test locally** (optional but recommended)
   ```bash
   conda build recipes/wasp2
   ```

5. **Submit PR**
   ```bash
   git checkout -b add-wasp2
   git add recipes/wasp2
   git commit -m "Add wasp2 recipe"
   git push origin add-wasp2
   gh pr create --repo bioconda/bioconda-recipes
   ```

## References

- [Bioconda Contributor Guide](https://bioconda.github.io/contributor/workflow.html)
- [Recipe Specification](https://bioconda.github.io/contributor/recipe-specification.html)
- [Example Recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes)

## Notes

- Windows builds are skipped (`skip: true  # [win]`) due to rust-htslib C dependencies
- The recipe requires both Rust and C compilers for the maturin build
- htslib is listed in both host and run requirements
