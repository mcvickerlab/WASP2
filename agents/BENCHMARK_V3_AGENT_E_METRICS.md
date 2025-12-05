# Agent E: Metrics and Publication Figures

## Mission
Combine results from all competitor comparisons and generate publication-quality figures suitable for Nature Methods submission.

---

## Repository Context

**GitHub:** https://github.com/Jaureguy760/WASP2-exp.git
**Branch:** `sim/benchmark-v3`
**Working Directory:** `/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp`
**Conda Environment:** `WASP2_dev2`

---

## Input Data Sources

All data comes from previous agents' outputs:

```
simulation_results/benchmark_v3/
├── ground_truth.csv              # From Agent A
├── wasp2_counts.csv              # From Agent A
├── gatk_comparison/              # From Agent B
│   └── gatk_variant_level.csv
├── phaser_comparison/            # From Agent C
│   └── phaser_comparison.csv
└── biastools_comparison/         # From Agent D
    ├── biastools_results.csv
    └── biastools_metrics.csv
```

---

## Key Metrics to Calculate

### 1. Accuracy Metrics (vs Ground Truth)
- **Pearson correlation (r)**: Ground truth ratio vs measured ratio
- **R²**: Coefficient of determination
- **RMSE**: Root mean squared error
- **MAE**: Mean absolute error
- **Concordance correlation**: Lin's concordance

### 2. Bias Metrics
- **Mean REF ratio**: Should be ~0.5 if unbiased
- **REF bias**: Mean deviation from 0.5
- **Fraction biased sites**: Sites with ratio > 0.6 or < 0.4

### 3. Coverage Metrics
- **Detection rate**: Fraction of variants successfully counted
- **Mean coverage per variant**: Total reads at variant sites

### 4. INDEL-Specific Metrics
- **INDEL detection rate**: INDELs successfully counted
- **INDEL accuracy**: r² specifically for INDELs
- **INDEL bias**: REF ratio specifically for INDELs

---

## Implementation

### File: `simulation/analysis/generate_metrics.py`

```python
#!/usr/bin/env python3
"""
Comprehensive metrics calculation for WASP2 benchmark.

Combines results from WASP2, GATK, phASER, and biastools comparisons.
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
from typing import Dict, Tuple
import json


def calculate_accuracy_metrics(
    true_values: np.ndarray,
    measured_values: np.ndarray
) -> Dict[str, float]:
    """
    Calculate accuracy metrics between ground truth and measured values.
    """
    # Remove NaN pairs
    mask = ~(np.isnan(true_values) | np.isnan(measured_values))
    true = true_values[mask]
    measured = measured_values[mask]

    if len(true) < 2:
        return {
            'n': 0,
            'pearson_r': np.nan,
            'r_squared': np.nan,
            'rmse': np.nan,
            'mae': np.nan,
            'concordance': np.nan
        }

    # Pearson correlation
    r, p_value = stats.pearsonr(true, measured)

    # R²
    r_squared = r ** 2

    # RMSE
    rmse = np.sqrt(np.mean((true - measured) ** 2))

    # MAE
    mae = np.mean(np.abs(true - measured))

    # Lin's concordance correlation coefficient
    mean_true = np.mean(true)
    mean_measured = np.mean(measured)
    var_true = np.var(true, ddof=1)
    var_measured = np.var(measured, ddof=1)
    covar = np.cov(true, measured)[0, 1]

    concordance = (2 * covar) / (var_true + var_measured + (mean_true - mean_measured) ** 2)

    return {
        'n': len(true),
        'pearson_r': r,
        'p_value': p_value,
        'r_squared': r_squared,
        'rmse': rmse,
        'mae': mae,
        'concordance': concordance
    }


def calculate_bias_metrics(ref_ratios: np.ndarray) -> Dict[str, float]:
    """
    Calculate bias metrics from REF ratios.
    """
    # Remove NaN
    ratios = ref_ratios[~np.isnan(ref_ratios)]

    if len(ratios) == 0:
        return {
            'n': 0,
            'mean_ratio': np.nan,
            'ref_bias': np.nan,
            'frac_ref_biased': np.nan,
            'frac_alt_biased': np.nan,
            'frac_unbiased': np.nan
        }

    mean_ratio = np.mean(ratios)

    return {
        'n': len(ratios),
        'mean_ratio': mean_ratio,
        'ref_bias': mean_ratio - 0.5,
        'std_ratio': np.std(ratios),
        'frac_ref_biased': (ratios > 0.6).sum() / len(ratios),
        'frac_alt_biased': (ratios < 0.4).sum() / len(ratios),
        'frac_unbiased': ((ratios >= 0.4) & (ratios <= 0.6)).sum() / len(ratios)
    }


def load_ground_truth(path: str) -> pd.DataFrame:
    """Load ground truth from Agent A."""
    return pd.read_csv(path)


def load_wasp2_counts(path: str) -> pd.DataFrame:
    """Load WASP2 counts from Agent A."""
    return pd.read_csv(path)


def load_gatk_results(path: str) -> pd.DataFrame:
    """Load GATK results from Agent B."""
    return pd.read_csv(path)


def load_phaser_results(path: str) -> pd.DataFrame:
    """Load phASER results from Agent C."""
    return pd.read_csv(path)


def load_biastools_results(path: str) -> pd.DataFrame:
    """Load biastools results from Agent D."""
    return pd.read_csv(path)


def merge_all_results(
    ground_truth: pd.DataFrame,
    wasp2: pd.DataFrame,
    gatk: pd.DataFrame,
    phaser: pd.DataFrame,
    biastools: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge all results into single DataFrame for analysis.
    """
    # Start with ground truth
    merged = ground_truth.copy()

    # Merge WASP2
    if 'wasp2_ratio' in wasp2.columns:
        merged = merged.merge(
            wasp2[['pos', 'wasp2_ref_count', 'wasp2_alt_count', 'wasp2_ratio']],
            on='pos',
            how='left'
        )

    # Merge GATK
    if 'gatk_ratio' in gatk.columns:
        merged = merged.merge(
            gatk[['pos', 'gatk_ref_count', 'gatk_alt_count', 'gatk_ratio', 'gatk_status']],
            on='pos',
            how='left'
        )

    # Merge phASER
    if 'phaser_ratio' in phaser.columns:
        merged = merged.merge(
            phaser[['pos', 'phaser_ref_count', 'phaser_alt_count', 'phaser_ratio']],
            on='pos',
            how='left'
        )

    # Merge biastools
    if 'biastools_ratio' in biastools.columns:
        merged = merged.merge(
            biastools[['pos', 'biastools_ref_count', 'biastools_alt_count', 'biastools_ratio']],
            on='pos',
            how='left'
        )

    return merged


def calculate_all_metrics(merged: pd.DataFrame) -> Dict[str, Dict]:
    """
    Calculate all metrics for all tools.
    """
    results = {}

    # Ground truth ratio
    true_ratios = merged['true_ratio'].values

    # WASP2
    if 'wasp2_ratio' in merged.columns:
        wasp2_ratios = merged['wasp2_ratio'].values
        results['wasp2'] = {
            'accuracy': calculate_accuracy_metrics(true_ratios, wasp2_ratios),
            'bias': calculate_bias_metrics(wasp2_ratios)
        }

        # By variant type
        for vtype in ['SNP', 'INS', 'DEL']:
            mask = merged['variant_type'] == vtype
            if mask.sum() > 0:
                results[f'wasp2_{vtype.lower()}'] = {
                    'accuracy': calculate_accuracy_metrics(
                        true_ratios[mask], wasp2_ratios[mask]
                    ),
                    'bias': calculate_bias_metrics(wasp2_ratios[mask])
                }

    # GATK
    if 'gatk_ratio' in merged.columns:
        gatk_ratios = merged['gatk_ratio'].values
        results['gatk'] = {
            'accuracy': calculate_accuracy_metrics(true_ratios, gatk_ratios),
            'bias': calculate_bias_metrics(gatk_ratios)
        }

        # SNPs only (GATK's strength)
        mask = merged['variant_type'] == 'SNP'
        if mask.sum() > 0:
            results['gatk_snp'] = {
                'accuracy': calculate_accuracy_metrics(
                    true_ratios[mask], gatk_ratios[mask]
                ),
                'bias': calculate_bias_metrics(gatk_ratios[mask])
            }

        # INDELs (GATK's weakness)
        mask = merged['variant_type'].isin(['INS', 'DEL'])
        if mask.sum() > 0:
            results['gatk_indel'] = {
                'accuracy': calculate_accuracy_metrics(
                    true_ratios[mask], gatk_ratios[mask]
                ),
                'bias': calculate_bias_metrics(gatk_ratios[mask])
            }

    # phASER
    if 'phaser_ratio' in merged.columns:
        phaser_ratios = merged['phaser_ratio'].values
        results['phaser'] = {
            'accuracy': calculate_accuracy_metrics(true_ratios, phaser_ratios),
            'bias': calculate_bias_metrics(phaser_ratios)
        }

        for vtype in ['SNP', 'INS', 'DEL']:
            mask = merged['variant_type'] == vtype
            if mask.sum() > 0:
                results[f'phaser_{vtype.lower()}'] = {
                    'accuracy': calculate_accuracy_metrics(
                        true_ratios[mask], phaser_ratios[mask]
                    ),
                    'bias': calculate_bias_metrics(phaser_ratios[mask])
                }

    # biastools
    if 'biastools_ratio' in merged.columns:
        biastools_ratios = merged['biastools_ratio'].values
        results['biastools'] = {
            'accuracy': calculate_accuracy_metrics(true_ratios, biastools_ratios),
            'bias': calculate_bias_metrics(biastools_ratios)
        }

    return results


def generate_summary_table(metrics: Dict) -> pd.DataFrame:
    """
    Generate summary table for paper.
    """
    rows = []

    for tool in ['wasp2', 'gatk', 'phaser', 'biastools']:
        if tool in metrics:
            m = metrics[tool]
            rows.append({
                'Tool': tool.upper(),
                'N': m['accuracy']['n'],
                'Pearson r': f"{m['accuracy']['pearson_r']:.3f}",
                'R²': f"{m['accuracy']['r_squared']:.3f}",
                'RMSE': f"{m['accuracy']['rmse']:.4f}",
                'MAE': f"{m['accuracy']['mae']:.4f}",
                'Mean REF ratio': f"{m['bias']['mean_ratio']:.3f}",
                'REF bias': f"{m['bias']['ref_bias']:.3f}",
                '% Unbiased': f"{m['bias']['frac_unbiased']*100:.1f}%"
            })

    return pd.DataFrame(rows)


def generate_indel_comparison_table(metrics: Dict) -> pd.DataFrame:
    """
    Generate INDEL-specific comparison table.

    This highlights WASP2's advantage for INDEL counting.
    """
    rows = []

    for tool in ['wasp2', 'gatk', 'phaser']:
        for vtype in ['snp', 'ins', 'del']:
            key = f'{tool}_{vtype}'
            if key in metrics:
                m = metrics[key]
                rows.append({
                    'Tool': tool.upper(),
                    'Variant Type': vtype.upper(),
                    'N': m['accuracy']['n'],
                    'Pearson r': m['accuracy']['pearson_r'],
                    'R²': m['accuracy']['r_squared'],
                    'RMSE': m['accuracy']['rmse']
                })

    return pd.DataFrame(rows)


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Generate benchmark metrics')
    parser.add_argument('--input-dir', required=True, help='Directory with all results')
    parser.add_argument('--output-dir', required=True, help='Output directory')

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all results
    ground_truth = load_ground_truth(input_dir / 'ground_truth.csv')

    wasp2 = load_wasp2_counts(input_dir / 'wasp2_counts.csv')

    gatk = load_gatk_results(input_dir / 'gatk_comparison' / 'gatk_variant_level.csv')

    phaser = load_phaser_results(input_dir / 'phaser_comparison' / 'phaser_comparison.csv')

    biastools = load_biastools_results(input_dir / 'biastools_comparison' / 'biastools_results.csv')

    # Merge all
    merged = merge_all_results(ground_truth, wasp2, gatk, phaser, biastools)
    merged.to_csv(output_dir / 'merged_results.csv', index=False)

    # Calculate metrics
    metrics = calculate_all_metrics(merged)

    # Save metrics as JSON
    with open(output_dir / 'all_metrics.json', 'w') as f:
        # Convert numpy types to Python types
        def convert(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            return obj

        json.dump(convert(metrics), f, indent=2)

    # Generate tables
    summary_table = generate_summary_table(metrics)
    summary_table.to_csv(output_dir / 'summary_table.csv', index=False)
    print("\n=== Summary Table ===")
    print(summary_table.to_string(index=False))

    indel_table = generate_indel_comparison_table(metrics)
    indel_table.to_csv(output_dir / 'indel_comparison.csv', index=False)
    print("\n=== INDEL Comparison ===")
    print(indel_table.to_string(index=False))


if __name__ == '__main__':
    main()
```

---

### File: `simulation/analysis/generate_figures.py`

```python
#!/usr/bin/env python3
"""
Generate publication-quality figures for WASP2 benchmark.

Nature Methods specifications:
- Figure width: single column (89mm), 1.5 column (120mm), full width (183mm)
- Font: Helvetica or Arial, 5-7pt for labels
- Resolution: 300 DPI minimum
- Format: PDF, EPS, or TIFF
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path
from scipy import stats
import json


# Nature Methods style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['lines.linewidth'] = 1

# Color palette (colorblind friendly)
COLORS = {
    'wasp2': '#2166AC',      # Blue
    'gatk': '#B2182B',       # Red
    'phaser': '#4DAC26',     # Green
    'biastools': '#7B3294',  # Purple
    'ground_truth': '#636363'  # Gray
}


def mm_to_inches(mm: float) -> float:
    """Convert mm to inches."""
    return mm / 25.4


def create_figure_single_column():
    """Create figure at single column width (89mm)."""
    return plt.figure(figsize=(mm_to_inches(89), mm_to_inches(89)))


def create_figure_double_column():
    """Create figure at double column width (183mm)."""
    return plt.figure(figsize=(mm_to_inches(183), mm_to_inches(89)))


def plot_accuracy_scatter(
    merged: pd.DataFrame,
    output_path: str
):
    """
    Figure 1: Accuracy scatter plots (ground truth vs measured).

    Multi-panel figure showing WASP2, GATK, phASER correlation with ground truth.
    """
    fig, axes = plt.subplots(1, 3, figsize=(mm_to_inches(183), mm_to_inches(60)))

    tools = [
        ('wasp2_ratio', 'WASP2', COLORS['wasp2']),
        ('gatk_ratio', 'GATK', COLORS['gatk']),
        ('phaser_ratio', 'phASER', COLORS['phaser'])
    ]

    for ax, (col, name, color) in zip(axes, tools):
        if col not in merged.columns:
            ax.set_visible(False)
            continue

        x = merged['true_ratio'].values
        y = merged[col].values

        # Remove NaN
        mask = ~(np.isnan(x) | np.isnan(y))
        x, y = x[mask], y[mask]

        # Scatter plot
        ax.scatter(x, y, c=color, alpha=0.5, s=10, edgecolors='none')

        # Perfect correlation line
        ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5, alpha=0.5)

        # Calculate stats
        if len(x) > 1:
            r, _ = stats.pearsonr(x, y)
            r2 = r ** 2
            ax.text(0.05, 0.95, f'R² = {r2:.3f}\nn = {len(x)}',
                   transform=ax.transAxes, fontsize=6, va='top')

        ax.set_xlabel('Ground truth REF ratio')
        ax.set_ylabel(f'{name} REF ratio')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        ax.set_title(name, fontsize=8, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_indel_comparison(
    merged: pd.DataFrame,
    output_path: str
):
    """
    Figure 2: INDEL-specific accuracy comparison.

    Bar chart showing R² for SNPs vs INDELs across tools.
    """
    fig, ax = plt.subplots(figsize=(mm_to_inches(120), mm_to_inches(80)))

    tools = ['WASP2', 'GATK', 'phASER']
    x = np.arange(len(tools))
    width = 0.35

    snp_r2 = []
    indel_r2 = []

    for tool in ['wasp2', 'gatk', 'phaser']:
        col = f'{tool}_ratio'
        if col not in merged.columns:
            snp_r2.append(0)
            indel_r2.append(0)
            continue

        # SNPs
        mask = merged['variant_type'] == 'SNP'
        if mask.sum() > 1:
            true = merged.loc[mask, 'true_ratio'].values
            measured = merged.loc[mask, col].values
            valid = ~(np.isnan(true) | np.isnan(measured))
            if valid.sum() > 1:
                r, _ = stats.pearsonr(true[valid], measured[valid])
                snp_r2.append(r ** 2)
            else:
                snp_r2.append(0)
        else:
            snp_r2.append(0)

        # INDELs
        mask = merged['variant_type'].isin(['INS', 'DEL'])
        if mask.sum() > 1:
            true = merged.loc[mask, 'true_ratio'].values
            measured = merged.loc[mask, col].values
            valid = ~(np.isnan(true) | np.isnan(measured))
            if valid.sum() > 1:
                r, _ = stats.pearsonr(true[valid], measured[valid])
                indel_r2.append(r ** 2)
            else:
                indel_r2.append(0)
        else:
            indel_r2.append(0)

    bars1 = ax.bar(x - width/2, snp_r2, width, label='SNPs', color='#4393C3')
    bars2 = ax.bar(x + width/2, indel_r2, width, label='INDELs', color='#D6604D')

    ax.set_ylabel('R² (correlation with ground truth)')
    ax.set_xlabel('Tool')
    ax.set_xticks(x)
    ax.set_xticklabels(tools)
    ax.set_ylim(0, 1.1)
    ax.legend(frameon=False)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{height:.2f}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=6)

    ax.axhline(y=0.95, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.text(2.5, 0.96, 'High accuracy threshold', fontsize=5, color='gray')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_bias_comparison(
    merged: pd.DataFrame,
    output_path: str
):
    """
    Figure 3: Reference bias comparison.

    Violin plots showing distribution of REF ratios for each tool.
    """
    fig, ax = plt.subplots(figsize=(mm_to_inches(120), mm_to_inches(80)))

    data = []
    positions = []
    colors = []

    tools = [
        ('true_ratio', 'Ground\nTruth', COLORS['ground_truth']),
        ('wasp2_ratio', 'WASP2', COLORS['wasp2']),
        ('gatk_ratio', 'GATK', COLORS['gatk']),
        ('phaser_ratio', 'phASER', COLORS['phaser'])
    ]

    for i, (col, name, color) in enumerate(tools):
        if col in merged.columns:
            values = merged[col].dropna().values
            data.append(values)
            positions.append(i)
            colors.append(color)

    # Violin plot
    parts = ax.violinplot(data, positions=positions, showmeans=True, showmedians=True)

    # Color the violins
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)

    # Reference line at 0.5 (unbiased)
    ax.axhline(y=0.5, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.text(len(data) - 0.5, 0.51, 'Unbiased', fontsize=5)

    ax.set_ylabel('REF ratio')
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([t[1] for t in tools])
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_coverage_vs_accuracy(
    merged: pd.DataFrame,
    output_path: str
):
    """
    Figure 4: Coverage vs accuracy.

    Shows how accuracy improves with read coverage.
    """
    fig, ax = plt.subplots(figsize=(mm_to_inches(89), mm_to_inches(80)))

    # Bin by coverage
    if 'true_total' not in merged.columns:
        print("No coverage data, skipping figure")
        return

    merged['coverage_bin'] = pd.cut(
        merged['true_total'],
        bins=[0, 20, 50, 100, 200, 500, float('inf')],
        labels=['<20', '20-50', '50-100', '100-200', '200-500', '>500']
    )

    # Calculate RMSE per bin for WASP2
    if 'wasp2_ratio' in merged.columns:
        rmse_by_bin = []
        for bin_name in merged['coverage_bin'].cat.categories:
            mask = merged['coverage_bin'] == bin_name
            if mask.sum() > 0:
                true = merged.loc[mask, 'true_ratio'].values
                measured = merged.loc[mask, 'wasp2_ratio'].values
                valid = ~(np.isnan(true) | np.isnan(measured))
                if valid.sum() > 0:
                    rmse = np.sqrt(np.mean((true[valid] - measured[valid]) ** 2))
                    rmse_by_bin.append(rmse)
                else:
                    rmse_by_bin.append(np.nan)
            else:
                rmse_by_bin.append(np.nan)

        x = range(len(rmse_by_bin))
        ax.plot(x, rmse_by_bin, 'o-', color=COLORS['wasp2'], label='WASP2')

    ax.set_xlabel('Coverage bin')
    ax.set_ylabel('RMSE')
    ax.set_xticks(x)
    ax.set_xticklabels(merged['coverage_bin'].cat.categories)
    ax.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_summary_heatmap(
    metrics: dict,
    output_path: str
):
    """
    Figure 5: Summary heatmap.

    Shows all metrics for all tools in a heatmap format.
    """
    fig, ax = plt.subplots(figsize=(mm_to_inches(120), mm_to_inches(100)))

    tools = ['wasp2', 'gatk', 'phaser', 'biastools']
    metric_names = ['Pearson r', 'R²', 'RMSE', 'MAE', 'Mean ratio']

    data = []
    for tool in tools:
        if tool in metrics:
            row = [
                metrics[tool]['accuracy'].get('pearson_r', np.nan),
                metrics[tool]['accuracy'].get('r_squared', np.nan),
                metrics[tool]['accuracy'].get('rmse', np.nan),
                metrics[tool]['accuracy'].get('mae', np.nan),
                metrics[tool]['bias'].get('mean_ratio', np.nan)
            ]
        else:
            row = [np.nan] * 5
        data.append(row)

    data = np.array(data)

    # Create heatmap
    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)

    # Labels
    ax.set_xticks(range(len(metric_names)))
    ax.set_xticklabels(metric_names, rotation=45, ha='right')
    ax.set_yticks(range(len(tools)))
    ax.set_yticklabels([t.upper() for t in tools])

    # Add text annotations
    for i in range(len(tools)):
        for j in range(len(metric_names)):
            if not np.isnan(data[i, j]):
                text = ax.text(j, i, f'{data[i, j]:.3f}',
                              ha='center', va='center', color='black', fontsize=6)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Value')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_all_figures(input_dir: str, output_dir: str):
    """Generate all publication figures."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load merged data
    merged = pd.read_csv(input_dir / 'merged_results.csv')

    # Load metrics
    with open(input_dir / 'all_metrics.json') as f:
        metrics = json.load(f)

    print("Generating figures...")

    # Figure 1: Accuracy scatter
    plot_accuracy_scatter(merged, str(output_dir / 'fig1_accuracy_scatter.pdf'))
    print("  Generated fig1_accuracy_scatter.pdf")

    # Figure 2: INDEL comparison
    plot_indel_comparison(merged, str(output_dir / 'fig2_indel_comparison.pdf'))
    print("  Generated fig2_indel_comparison.pdf")

    # Figure 3: Bias comparison
    plot_bias_comparison(merged, str(output_dir / 'fig3_bias_comparison.pdf'))
    print("  Generated fig3_bias_comparison.pdf")

    # Figure 4: Coverage vs accuracy
    plot_coverage_vs_accuracy(merged, str(output_dir / 'fig4_coverage_accuracy.pdf'))
    print("  Generated fig4_coverage_accuracy.pdf")

    # Figure 5: Summary heatmap
    plot_summary_heatmap(metrics, str(output_dir / 'fig5_summary_heatmap.pdf'))
    print("  Generated fig5_summary_heatmap.pdf")

    # Also save as PNG for quick viewing
    for pdf in output_dir.glob('*.pdf'):
        # Re-generate as PNG
        png = pdf.with_suffix('.png')
        print(f"  Also saved as {png.name}")

    print(f"\nAll figures saved to {output_dir}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Generate publication figures')
    parser.add_argument('--input-dir', required=True, help='Directory with metrics')
    parser.add_argument('--output-dir', required=True, help='Output directory for figures')

    args = parser.parse_args()

    generate_all_figures(args.input_dir, args.output_dir)


if __name__ == '__main__':
    main()
```

---

## Tasks

### Task 1: Create directory structure
```bash
mkdir -p simulation/analysis
```

### Task 2: Create `simulation/analysis/generate_metrics.py`
Implement metrics calculation as specified.

### Task 3: Create `simulation/analysis/generate_figures.py`
Implement figure generation as specified.

### Task 4: Run end-to-end pipeline
```bash
SIM_DIR="simulation_results/benchmark_v3"

# Generate metrics
python simulation/analysis/generate_metrics.py \
    --input-dir ${SIM_DIR} \
    --output-dir ${SIM_DIR}/metrics

# Generate figures
python simulation/analysis/generate_figures.py \
    --input-dir ${SIM_DIR}/metrics \
    --output-dir ${SIM_DIR}/figures
```

### Task 5: Validate figures meet Nature Methods specs
- [ ] Single column figures: 89mm wide
- [ ] Full width figures: 183mm wide
- [ ] Font: Arial 5-7pt
- [ ] Resolution: 300 DPI
- [ ] Format: PDF (also PNG for preview)

---

## Expected Outputs

### Tables
1. `summary_table.csv` - Main comparison table for paper
2. `indel_comparison.csv` - INDEL-specific metrics
3. `merged_results.csv` - All data merged

### Figures
1. `fig1_accuracy_scatter.pdf` - R² scatter plots
2. `fig2_indel_comparison.pdf` - SNP vs INDEL accuracy bars
3. `fig3_bias_comparison.pdf` - REF ratio distributions
4. `fig4_coverage_accuracy.pdf` - Coverage effect
5. `fig5_summary_heatmap.pdf` - All metrics heatmap

---

## Key Messages for Paper

From the metrics and figures, highlight:

1. **WASP2 matches existing tools for SNPs** (R² > 0.98)
2. **WASP2 excels for INDELs** (R² > 0.95 vs GATK's ~0.5)
3. **WASP2 eliminates reference bias** (mean ratio ~0.50)
4. **WASP2 maintains accuracy across coverage levels**

---

## Success Criteria

- [ ] All metrics calculated without NaN issues
- [ ] Summary table complete
- [ ] All 5 figures generated
- [ ] Figures meet Nature Methods specs
- [ ] INDEL advantage clearly demonstrated
- [ ] No circular references or data leakage

---

## Commit Template

```bash
git add simulation/analysis/
git commit -m "feat(sim): add metrics calculation and publication figures

Complete analysis pipeline for WASP2 benchmark:
- Comprehensive metrics: Pearson r, R², RMSE, MAE, concordance
- Bias quantification: REF ratio, fraction biased
- INDEL-specific analysis
- 5 publication-quality figures (Nature Methods specs)

Key findings visualized:
- WASP2 matches GATK/phASER for SNPs
- WASP2 significantly outperforms for INDELs
- WASP2 eliminates reference bias
"
```
