# Agent 4: Publication Metrics and Figures

## Mission
Compute publication-quality metrics and generate figures for Nature Methods submission. This agent runs **after** Agents 1-3 complete and produce data.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/metrics`
**Parent Branch:** `ropc-indels`

**Working Directory:**
```
/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp
```

**Conda Environment:** `WASP2_dev2`

---

## Dependencies

**WAIT for these to complete before starting:**

| Agent | Branch | Output | Status |
|-------|--------|--------|--------|
| 1 | sim/comprehensive | `simulation_results/comprehensive_*/simulation_results.csv` | Required |
| 2 | sim/paired-end | `simulation_results/paired_end_*/simulation_results.csv` | Required |
| 3 | sim/gatk-compare | `comparison_results/merged_counts.csv` | Required |

---

## Inputs Required

### From Agent 1 (Comprehensive):
```
simulation_results/comprehensive_TIMESTAMP/
├── simulation_results.csv    # 810 tests with ground truth
└── reference.fa              # For variant context
```

**CSV Columns:**
```
chrom, pos, variant_type, coverage, replicate, true_ratio,
ref_count, alt_count, total_reads, observed_ratio, error, error_pct, status
```

### From Agent 2 (Paired-End):
```
simulation_results/paired_end_TIMESTAMP/
└── simulation_results.csv    # Same format as above
```

### From Agent 3 (GATK Compare):
```
comparison_results/
├── merged_counts.csv         # WASP2 + GATK + truth counts
├── gatk_metrics.json         # GATK accuracy metrics
└── wasp2_metrics.json        # WASP2 accuracy metrics
```

**merged_counts.csv Columns:**
```
chrom, pos, truth_ref, truth_alt, truth_total, truth_ratio, variant_type,
gatk_ref, gatk_alt, gatk_total, wasp2_ref, wasp2_alt, wasp2_total
```

---

## Publication Figure Specifications

### Nature Methods Requirements:
- **Resolution:** 300 DPI minimum
- **Formats:** PDF (vector) + PNG (raster)
- **Font:** Arial or Helvetica, 6-8pt minimum
- **Colors:** Colorblind-friendly palette
- **Width:** Single column (86mm) or double column (178mm)

### Recommended Color Palette:
```python
COLORS = {
    'SNP': '#2ecc71',      # Green
    'INS': '#3498db',      # Blue
    'DEL': '#e74c3c',      # Red
    'WASP2': '#3498db',    # Blue
    'GATK': '#e74c3c',     # Red
    'truth': '#2c3e50',    # Dark gray
}
```

---

## Figures to Generate

### Figure 1: Ground Truth Correlation (Main Result)
**Purpose:** Show WASP2 accurately recovers allelic ratios

```python
def generate_figure_1(results_df: pd.DataFrame, output_dir: Path):
    """
    Scatter plot: Observed vs Expected allelic ratio
    - X-axis: True ratio (ground truth)
    - Y-axis: Observed ratio (WASP2)
    - Color by variant type (SNP/INS/DEL)
    - Include Pearson r and regression line
    - Size: 6x6 inches (single column, square)
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    # Filter infinite values
    valid = results_df[results_df['observed_ratio'] != np.inf].copy()

    # Plot by variant type
    for vtype, color in [('SNP', '#2ecc71'), ('INS', '#3498db'), ('DEL', '#e74c3c')]:
        subset = valid[valid['variant_type'] == vtype]
        ax.scatter(
            subset['true_ratio'],
            subset['observed_ratio'],
            c=color,
            label=f'{vtype} (n={len(subset)})',
            alpha=0.6,
            s=30,
            edgecolors='white',
            linewidth=0.5
        )

    # Perfect correlation line
    max_val = max(valid['true_ratio'].max(), valid['observed_ratio'].max()) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1, label='y = x')

    # Compute statistics
    r, p = pearsonr(valid['observed_ratio'], valid['true_ratio'])

    # Labels and formatting
    ax.set_xlabel('Expected Allelic Ratio (REF/ALT)', fontsize=10)
    ax.set_ylabel('Observed Allelic Ratio (REF/ALT)', fontsize=10)
    ax.set_title(f'WASP2 Allelic Ratio Accuracy\nPearson r = {r:.4f}', fontsize=12)
    ax.legend(loc='upper left', fontsize=8)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal')

    # Add grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(output_dir / 'figure1_correlation.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure1_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Figure 1 saved: r = {r:.4f}")
```

### Figure 2: Error Distribution by Variant Type
**Purpose:** Show consistent accuracy across variant types

```python
def generate_figure_2(results_df: pd.DataFrame, output_dir: Path):
    """
    Box plot: Error percentage by variant type
    - X-axis: Variant type (SNP, INS, DEL)
    - Y-axis: Error (%)
    - Show individual points
    - Include 10% threshold line
    - Size: 8x5 inches
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    order = ['SNP', 'INS', 'DEL']
    colors = ['#2ecc71', '#3498db', '#e74c3c']

    # Box plot
    bp = ax.boxplot(
        [results_df[results_df['variant_type'] == vt]['error_pct'] for vt in order],
        labels=order,
        patch_artist=True,
        widths=0.6
    )

    # Color boxes
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add individual points (jittered)
    for i, vtype in enumerate(order):
        subset = results_df[results_df['variant_type'] == vtype]
        x = np.random.normal(i + 1, 0.08, size=len(subset))
        ax.scatter(x, subset['error_pct'], alpha=0.4, s=15, c='black')

    # 10% threshold
    ax.axhline(y=10, color='red', linestyle='--', alpha=0.7, linewidth=1.5,
               label='10% threshold (PASS/FAIL)')

    # Formatting
    ax.set_xlabel('Variant Type', fontsize=10)
    ax.set_ylabel('Error (%)', fontsize=10)
    ax.set_title('WASP2 Error Distribution by Variant Type', fontsize=12)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_ylim(bottom=0)
    ax.grid(True, axis='y', alpha=0.3)

    # Add pass rates as text
    for i, vtype in enumerate(order):
        subset = results_df[results_df['variant_type'] == vtype]
        pass_rate = (subset['status'] == 'PASS').mean() * 100
        ax.text(i + 1, ax.get_ylim()[1] * 0.95, f'{pass_rate:.0f}%',
                ha='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_dir / 'figure2_error_by_type.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure2_error_by_type.png', dpi=300, bbox_inches='tight')
    plt.close()
```

### Figure 3: WASP2 vs GATK Comparison
**Purpose:** Show WASP2 is comparable or better than GATK

```python
def generate_figure_3(comparison_df: pd.DataFrame, output_dir: Path):
    """
    Side-by-side scatter: WASP2 vs GATK accuracy
    - Left panel: WASP2 vs truth
    - Right panel: GATK vs truth
    - Same scale for comparison
    - Size: 12x5 inches (double column)
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Compute ratios
    df = comparison_df.copy()
    df['truth_ratio'] = df['truth_ref'] / df['truth_total']
    df['wasp2_ratio'] = df['wasp2_ref'] / df['wasp2_total']
    df['gatk_ratio'] = df['gatk_ref'] / df['gatk_total']

    # Filter valid
    df = df.dropna()

    max_val = 1.1

    # WASP2 panel
    ax1 = axes[0]
    r_wasp2, _ = pearsonr(df['wasp2_ratio'], df['truth_ratio'])
    ax1.scatter(df['truth_ratio'], df['wasp2_ratio'], alpha=0.5, s=20, c='#3498db')
    ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    ax1.set_xlabel('True REF/(REF+ALT)', fontsize=10)
    ax1.set_ylabel('WASP2 REF/(REF+ALT)', fontsize=10)
    ax1.set_title(f'WASP2 (r = {r_wasp2:.4f})', fontsize=12)
    ax1.set_xlim(0, max_val)
    ax1.set_ylim(0, max_val)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)

    # GATK panel
    ax2 = axes[1]
    r_gatk, _ = pearsonr(df['gatk_ratio'], df['truth_ratio'])
    ax2.scatter(df['truth_ratio'], df['gatk_ratio'], alpha=0.5, s=20, c='#e74c3c')
    ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    ax2.set_xlabel('True REF/(REF+ALT)', fontsize=10)
    ax2.set_ylabel('GATK REF/(REF+ALT)', fontsize=10)
    ax2.set_title(f'GATK ASEReadCounter (r = {r_gatk:.4f})', fontsize=12)
    ax2.set_xlim(0, max_val)
    ax2.set_ylim(0, max_val)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)

    plt.suptitle('Allelic Ratio Accuracy: WASP2 vs GATK', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'figure3_wasp2_vs_gatk.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure3_wasp2_vs_gatk.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Figure 3 saved: WASP2 r={r_wasp2:.4f}, GATK r={r_gatk:.4f}")
```

### Figure 4: Coverage Effect
**Purpose:** Show consistent accuracy across coverage levels

```python
def generate_figure_4(results_df: pd.DataFrame, output_dir: Path):
    """
    Line plot: Error vs coverage level
    - X-axis: Coverage (20x, 50x, 100x)
    - Y-axis: Mean error (%)
    - Lines for each variant type
    - Error bars for std dev
    - Size: 7x5 inches
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    coverages = sorted(results_df['coverage'].unique())
    colors = {'SNP': '#2ecc71', 'INS': '#3498db', 'DEL': '#e74c3c'}

    for vtype in ['SNP', 'INS', 'DEL']:
        means = []
        stds = []
        for cov in coverages:
            subset = results_df[(results_df['variant_type'] == vtype) &
                                (results_df['coverage'] == cov)]
            means.append(subset['error_pct'].mean())
            stds.append(subset['error_pct'].std())

        ax.errorbar(coverages, means, yerr=stds, marker='o', label=vtype,
                    color=colors[vtype], linewidth=2, markersize=8, capsize=5)

    ax.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='10% threshold')
    ax.set_xlabel('Coverage (x)', fontsize=10)
    ax.set_ylabel('Mean Error (%)', fontsize=10)
    ax.set_title('WASP2 Accuracy vs Coverage', fontsize=12)
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xticks(coverages)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(output_dir / 'figure4_coverage_effect.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure4_coverage_effect.png', dpi=300, bbox_inches='tight')
    plt.close()
```

---

## Metrics to Compute

### Summary Table for Paper:

```python
def compute_summary_table(results_df: pd.DataFrame) -> pd.DataFrame:
    """Generate summary statistics table for paper."""

    summary = []

    # Overall
    summary.append({
        'Category': 'Overall',
        'N': len(results_df),
        'Pass Rate (%)': (results_df['status'] == 'PASS').mean() * 100,
        'Mean Error (%)': results_df['error_pct'].mean(),
        'Median Error (%)': results_df['error_pct'].median(),
        'Max Error (%)': results_df['error_pct'].max(),
    })

    # By variant type
    for vtype in ['SNP', 'INS', 'DEL']:
        subset = results_df[results_df['variant_type'] == vtype]
        summary.append({
            'Category': vtype,
            'N': len(subset),
            'Pass Rate (%)': (subset['status'] == 'PASS').mean() * 100,
            'Mean Error (%)': subset['error_pct'].mean(),
            'Median Error (%)': subset['error_pct'].median(),
            'Max Error (%)': subset['error_pct'].max(),
        })

    return pd.DataFrame(summary)
```

---

## Full Implementation

### File: `simulation/generate_publication_figures.py`

```python
#!/usr/bin/env python3
"""
Generate Publication Figures for WASP2 Simulation Benchmarks

Requires outputs from:
- Agent 1 (sim/comprehensive)
- Agent 2 (sim/paired-end)
- Agent 3 (sim/gatk-compare)

Usage:
    python generate_publication_figures.py \
        --comprehensive simulation_results/comprehensive_*/simulation_results.csv \
        --paired-end simulation_results/paired_end_*/simulation_results.csv \
        --comparison comparison_results/merged_counts.csv \
        --output figures/
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 8

# [Include all figure generation functions from above]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comprehensive', required=True)
    parser.add_argument('--paired-end', help='Optional paired-end results')
    parser.add_argument('--comparison', help='Optional GATK comparison')
    parser.add_argument('--output', required=True)

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("Generating Publication Figures")
    print("="*60)

    # Load comprehensive results
    print("\nLoading comprehensive results...")
    comp_df = pd.read_csv(args.comprehensive)
    print(f"  {len(comp_df)} tests")

    # Figure 1: Correlation
    print("\nGenerating Figure 1 (correlation)...")
    generate_figure_1(comp_df, output_dir)

    # Figure 2: Error by type
    print("Generating Figure 2 (error by type)...")
    generate_figure_2(comp_df, output_dir)

    # Figure 4: Coverage effect
    if 'coverage' in comp_df.columns and comp_df['coverage'].nunique() > 1:
        print("Generating Figure 4 (coverage effect)...")
        generate_figure_4(comp_df, output_dir)

    # Figure 3: GATK comparison (if available)
    if args.comparison and Path(args.comparison).exists():
        print("\nLoading GATK comparison...")
        compare_df = pd.read_csv(args.comparison)
        print("Generating Figure 3 (WASP2 vs GATK)...")
        generate_figure_3(compare_df, output_dir)

    # Summary table
    print("\nGenerating summary table...")
    summary = compute_summary_table(comp_df)
    summary.to_csv(output_dir / 'summary_table.csv', index=False)
    print(summary.to_string(index=False))

    # Compute publication metrics
    metrics = {
        'n_tests': len(comp_df),
        'pass_rate': (comp_df['status'] == 'PASS').mean() * 100,
        'mean_error': comp_df['error_pct'].mean(),
        'median_error': comp_df['error_pct'].median(),
    }

    # Correlation
    valid = comp_df[comp_df['observed_ratio'] != np.inf]
    if len(valid) > 0:
        metrics['pearson_r'], metrics['pearson_p'] = pearsonr(
            valid['observed_ratio'], valid['true_ratio']
        )

    with open(output_dir / 'publication_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"Figures saved to: {output_dir}")
    print(f"\nKey metrics:")
    print(f"  Pass rate: {metrics['pass_rate']:.1f}%")
    print(f"  Pearson r: {metrics.get('pearson_r', 'N/A'):.4f}")
    print(f"  Mean error: {metrics['mean_error']:.2f}%")


if __name__ == '__main__':
    main()
```

---

## Step-by-Step Execution

### Step 1: Wait for Dependencies
```bash
# Check if Agent 1 completed
ls simulation_results/comprehensive_*/simulation_results.csv

# Check if Agent 3 completed
ls comparison_results/merged_counts.csv
```

### Step 2: Setup
```bash
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate WASP2_dev2
cd /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp

git checkout sim/metrics
git pull origin sim/metrics

# Merge completed branches
git merge sim/comprehensive --no-edit
git merge sim/gatk-compare --no-edit
# git merge sim/paired-end --no-edit  # if available

mkdir -p simulation figures
```

### Step 3: Create and Run
```bash
# Create simulation/generate_publication_figures.py

# Run
COMP_CSV=$(ls -t simulation_results/comprehensive_*/simulation_results.csv | head -1)
COMPARE_CSV="comparison_results/merged_counts.csv"

python simulation/generate_publication_figures.py \
    --comprehensive ${COMP_CSV} \
    --comparison ${COMPARE_CSV} \
    --output figures/
```

### Step 4: Validate Outputs
```bash
ls figures/
# Should see:
# - figure1_correlation.pdf
# - figure1_correlation.png
# - figure2_error_by_type.pdf
# - figure2_error_by_type.png
# - figure3_wasp2_vs_gatk.pdf (if GATK comparison available)
# - figure3_wasp2_vs_gatk.png
# - figure4_coverage_effect.pdf
# - figure4_coverage_effect.png
# - summary_table.csv
# - publication_metrics.json
```

---

## Success Criteria

- [ ] Figure 1 (correlation): Pearson r > 0.95, clean scatter plot
- [ ] Figure 2 (error by type): All variant types < 10% mean error
- [ ] Figure 3 (WASP2 vs GATK): WASP2 r ≥ GATK r
- [ ] Figure 4 (coverage): Consistent across coverage levels
- [ ] Summary table generated
- [ ] All figures in PDF + PNG format at 300 DPI

---

## Commit Template

```bash
git add simulation/generate_publication_figures.py
git add figures/

git commit -m "feat: add publication figures and metrics

Figures generated:
- Figure 1: Ground truth correlation (r = X.XXXX)
- Figure 2: Error distribution by variant type
- Figure 3: WASP2 vs GATK comparison
- Figure 4: Coverage effect analysis

Summary:
- Total tests: 810
- Overall pass rate: XX%
- SNP: XX% pass, X.X% mean error
- INS: XX% pass, X.X% mean error
- DEL: XX% pass, X.X% mean error
"

git push origin sim/metrics
```

---

## Final Merge to ropc-indels

After all figures are validated:

```bash
git checkout ropc-indels
git merge sim/metrics --no-edit
git push origin ropc-indels
```

---

## Handoff

When complete:
1. All figures in `figures/` directory
2. Summary metrics in `publication_metrics.json`
3. Ready for paper writing
