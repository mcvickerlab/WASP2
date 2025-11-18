# WASP2 Development Container

This directory contains the configuration for GitHub Codespaces and VS Code Dev Containers.

## 🚀 Quick Start with GitHub Codespaces

### Option 1: Via GitHub Web Interface

1. Go to your repository: `https://github.com/Jaureguy760/WASP2-exp`
2. Click the green **"Code"** button
3. Select the **"Codespaces"** tab
4. Click **"Create codespace on [branch-name]"**
5. Wait 2-3 minutes for the environment to build
6. Start coding! Everything is pre-installed.

### Option 2: Via GitHub CLI

```bash
# Create a new codespace
gh codespace create --repo Jaureguy760/WASP2-exp

# List your codespaces
gh codespace list

# Connect to a codespace
gh codespace code
```

### Option 3: Via VS Code Desktop

1. Install the **GitHub Codespaces** extension in VS Code
2. Press `F1` → "Codespaces: Create New Codespace"
3. Select `Jaureguy760/WASP2-exp`
4. Select branch
5. Environment builds automatically

## 🛠️ What Gets Installed

The devcontainer automatically sets up:

### System Tools
- ✅ Conda/Mamba package manager
- ✅ Git & GitHub CLI
- ✅ Python 3.11

### Bioinformatics Tools
- ✅ bedtools
- ✅ bcftools
- ✅ samtools
- ✅ pysam
- ✅ pybedtools

### Python Packages
- ✅ All dependencies from `environment.yml`
- ✅ All dependencies from `requirements.txt`
- ✅ Development tools (pytest, mypy, black)
- ✅ WASP2 installed in editable mode (`pip install -e .`)

### VS Code Extensions
- ✅ Python language support
- ✅ Pylance (type checking)
- ✅ Mypy type checker
- ✅ Black formatter
- ✅ Jupyter notebooks
- ✅ TOML support
- ✅ GitHub Copilot (if you have access)

### Pre-configured Settings
- ✅ Auto-format on save (Black)
- ✅ Type checking enabled (mypy)
- ✅ Pytest test discovery
- ✅ 88-character line ruler
- ✅ Pre-commit hooks installed

## 📋 Using the Environment

### Run Tests
```bash
pytest tests/
pytest tests/regression/
```

### Type Checking
```bash
mypy src/
```

### Format Code
```bash
black src/
# Or just save a file - auto-format is enabled!
```

### Run WASP2 Commands
```bash
wasp2-count --help
wasp2-map --help
wasp2-analyze --help
```

### Run Full Pipeline
```bash
# Example with test data
bash scripts/run_full_pipeline_baseline.sh
```

## 🔧 Local Dev Container (VS Code)

If you want to run the dev container locally instead of in the cloud:

### Requirements
- Docker Desktop installed
- VS Code with "Dev Containers" extension

### Steps
1. Clone the repo: `git clone https://github.com/Jaureguy760/WASP2-exp.git`
2. Open in VS Code: `code WASP2-exp`
3. Press `F1` → "Dev Containers: Reopen in Container"
4. Wait for build (first time takes 5-10 minutes)
5. Container is ready!

## 🎯 Benefits

### For Contributors
- **Zero setup time** - no need to install conda, bioinformatics tools, etc.
- **Consistent environment** - everyone uses the exact same tools and versions
- **Pre-configured** - VS Code settings, extensions, and tools ready to go
- **Isolated** - doesn't mess with your local machine

### For Maintainers
- **Easier onboarding** - new contributors can start coding in minutes
- **Reproducible** - issues can be debugged in identical environments
- **Version controlled** - dev environment config is in git
- **Cloud-based** - work from anywhere, any device

## 📁 Files

- **`devcontainer.json`** - Main configuration file
  - Defines base image (miniconda)
  - Lists VS Code extensions to install
  - Configures editor settings
  - Runs setup script on creation

- **`setup.sh`** - Post-create setup script
  - Creates conda environment from `environment.yml`
  - Installs WASP2 in editable mode
  - Sets up pre-commit hooks
  - Verifies installation

- **`README.md`** - This file!

## 🔄 Updating the Environment

If dependencies change (new packages added to `environment.yml` or `requirements.txt`):

### In an existing Codespace:
```bash
# Update conda environment
conda env update -f environment.yml --prune

# Reinstall package
pip install -e .
```

### For a fresh start:
1. Delete the old codespace
2. Create a new one (will use updated config)

## 🐛 Troubleshooting

### Codespace won't start
- Check GitHub status: https://www.githubstatus.com/
- Try deleting and recreating the codespace
- Check repository permissions

### Missing packages
```bash
# Rebuild conda environment
conda env create -f environment.yml --force
```

### VS Code extensions not loading
- Reload window: `F1` → "Developer: Reload Window"
- Rebuild container: `F1` → "Codespaces: Rebuild Container"

### Pre-commit hooks failing
```bash
# Reinstall hooks
pre-commit uninstall
pre-commit install
```

## 💡 Tips

1. **Preserve your work**: Codespaces auto-save, but commit frequently!
2. **Port forwarding**: Automatically forwards ports if you run servers
3. **Dotfiles**: GitHub Codespaces can apply your personal dotfiles
4. **Secrets**: Use GitHub Secrets for API keys (Settings → Secrets → Codespaces)
5. **Multiple codespaces**: You can have multiple for different branches

## 🔗 Resources

- [GitHub Codespaces Docs](https://docs.github.com/en/codespaces)
- [VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/containers)
- [Dev Container Spec](https://containers.dev/)
- [WASP2 Documentation](https://wasp2.readthedocs.io)

## 💰 Costs

### GitHub Codespaces Pricing
- **Free tier**: 120 core-hours/month for personal accounts
- **Pro accounts**: More free hours included
- **After free tier**: ~$0.18/hour for 2-core machines

### Tips to save hours:
- Stop codespace when not using (auto-stops after 30 min idle)
- Delete unused codespaces
- Use VS Code desktop to connect (doesn't count as browser hours)

---

**Questions?** Check the [GitHub Codespaces docs](https://docs.github.com/en/codespaces) or open an issue!
