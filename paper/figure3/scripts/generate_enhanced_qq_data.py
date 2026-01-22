#!/usr/bin/env python3
"""
Enhanced QQ Data Generation for Method Comparison

Aggregates p-values from 137 RNA-seq samples across 24+ statistical methods
to generate publication-quality QQ plot data with:
- Genomic inflation factor (lambda) per method
- Method metadata (dispersion type, FDR method)
- Confidence bands for expected distribution

Data source:
- cvpc/experiments/stage1_global_dispersion/results/ (exp01-exp09, exp12-20)
- cvpc/experiments/stage1_cov_dispersion_bins/results/ (exp10, exp13-17)

Methods tested:
1. Binomial (naive, no dispersion)
2. Beta-binomial with global dispersion (MoM, MLE)
3. Beta-binomial with coverage-binned dispersion
4. Various FDR methods (BH, Discrete, Mid-p)

Author: WASP2 R&D Team
"""
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from scipy import stats
import warnings

# Suppress pandas SettingWithCopyWarning
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)

# Add paper directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config import REPO_ROOT

# Data paths
CVPC_ROOT = Path("/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc")
GLOBAL_DISP_DIR = CVPC_ROOT / "experiments/stage1_global_dispersion/results"
COV_BINS_DIR = CVPC_ROOT / "experiments/stage1_cov_dispersion_bins/results"

# Method metadata - describes each experiment
METHOD_METADATA = {
    'exp01': {'name': 'Binomial', 'dispersion': 'None', 'fdr': 'BH', 'estimator': 'naive'},
    'exp02': {'name': 'BB-WASP2-MLE', 'dispersion': 'Global', 'fdr': 'BH', 'estimator': 'WASP2 MLE'},
    'exp03': {'name': 'BB-Custom-MLE', 'dispersion': 'Global', 'fdr': 'BH', 'estimator': 'Custom MLE'},
    'exp04': {'name': 'BB-MoM', 'dispersion': 'Global', 'fdr': 'BH', 'estimator': 'MoM'},
    'exp05': {'name': 'BB-WASP2-Discrete', 'dispersion': 'Global', 'fdr': 'Discrete', 'estimator': 'WASP2 MLE'},
    'exp06': {'name': 'BB-Custom-Discrete', 'dispersion': 'Global', 'fdr': 'Discrete', 'estimator': 'Custom MLE'},
    'exp07': {'name': 'BB-MoM-Discrete', 'dispersion': 'Global', 'fdr': 'Discrete', 'estimator': 'MoM'},
    'exp08': {'name': 'BB-Custom-Midp', 'dispersion': 'Global', 'fdr': 'Mid-p', 'estimator': 'Custom MLE'},
    'exp09': {'name': 'BB-MoM-Midp', 'dispersion': 'Global', 'fdr': 'Mid-p', 'estimator': 'MoM'},
    'exp10': {'name': 'BB-CovBins', 'dispersion': 'Coverage-binned', 'fdr': 'Discrete', 'estimator': 'Binned'},
    'exp12': {'name': 'BB-Lowess', 'dispersion': 'Coverage-LOWESS', 'fdr': 'BH', 'estimator': 'LOWESS'},
    'exp13': {'name': 'BB-CovBins-BH', 'dispersion': 'Coverage-binned', 'fdr': 'BH', 'estimator': 'Binned MoM'},
    'exp15': {'name': 'BB-CovBins-MLE', 'dispersion': 'Coverage-binned', 'fdr': 'BH', 'estimator': 'Binned MLE'},
    'exp16': {'name': 'BB-CovBins-MLE-Disc', 'dispersion': 'Coverage-binned', 'fdr': 'Discrete', 'estimator': 'Binned MLE'},
    'exp17': {'name': 'BB-CovBins-MoM-Midp', 'dispersion': 'Coverage-binned', 'fdr': 'Mid-p', 'estimator': 'Binned MoM'},
    'exp20': {'name': 'BB-Linear', 'dispersion': 'Coverage-linear', 'fdr': 'BH', 'estimator': 'Linear'},
}


def compute_lambda_gc(pvalues: np.ndarray) -> float:
    """
    Compute genomic inflation factor (lambda_GC).

    Lambda measures deviation from expected uniform distribution:
    - lambda = 1.0: well-calibrated (ideal)
    - lambda > 1.0: inflation (too many small p-values)
    - lambda < 1.0: deflation (too many large p-values)

    Formula: lambda = median(chi^2_observed) / median(chi^2_expected)
             = median(-2*log(p)) / 0.455

    Args:
        pvalues: Array of p-values

    Returns:
        Genomic inflation factor
    """
    # Filter valid p-values
    pv = pvalues[(pvalues > 0) & (pvalues <= 1) & ~np.isnan(pvalues)]

    if len(pv) < 100:
        return np.nan

    # Convert to chi-squared statistics
    chi2_obs = stats.chi2.isf(pv, df=1)
    median_obs = np.median(chi2_obs)

    # Expected median for chi2(1) is 0.455
    median_exp = stats.chi2.ppf(0.5, df=1)  # = 0.4549...

    return median_obs / median_exp


def load_experiment_pvalues(exp_id: str, results_dir: Path, max_samples: int = None) -> np.ndarray:
    """
    Load and aggregate p-values from all samples for a given experiment.

    Args:
        exp_id: Experiment identifier (e.g., 'exp01', 'exp15')
        results_dir: Directory containing result TSV files
        max_samples: Maximum samples to load (for testing)

    Returns:
        Array of all p-values across samples
    """
    # Try multiple naming patterns
    # exp01_* for basic experiments
    # exp15_b10_* for coverage-binned experiments (b10 = 10 bins)
    # exp12_lo10_* for LOWESS experiments (lo10 = 10 iterations)
    patterns = [
        f"{exp_id}_*.tsv",
        f"{exp_id}_b10_*.tsv",  # Coverage bins with 10 bins
        f"{exp_id}_lo10_*.tsv",  # LOWESS with 10 iterations
        f"{exp_id}_l10_*.tsv",   # Linear with 10 bins
    ]

    files = []
    for pattern in patterns:
        found = sorted(results_dir.glob(pattern))
        # Filter out dispersion parameter files (they have 'dispersion_bins' in name)
        found = [f for f in found if 'dispersion_bins' not in f.name]
        if found:
            files = found
            break

    if not files:
        print(f"  WARNING: No files found for {exp_id} in {results_dir}")
        return np.array([])

    if max_samples:
        files = files[:max_samples]

    all_pvalues = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t', usecols=['pvalue'])
            pv = df['pvalue'].dropna().values
            all_pvalues.extend(pv)
        except Exception as e:
            print(f"  Warning: Error reading {f}: {e}")

    return np.array(all_pvalues)


def generate_qq_points(pvalues: np.ndarray, n_points: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate QQ plot points (expected vs observed p-values).

    Uses quantile thinning to reduce to n_points for plotting efficiency.

    Args:
        pvalues: Array of observed p-values
        n_points: Number of points to return

    Returns:
        (expected, observed) arrays in -log10 scale
    """
    pv = pvalues[(pvalues > 0) & (pvalues <= 1) & ~np.isnan(pvalues)]
    pv = np.sort(pv)

    n = len(pv)
    if n == 0:
        return np.array([]), np.array([])

    # Expected p-values under uniform distribution
    expected = np.arange(1, n + 1) / (n + 1)

    # Thin to n_points using quantiles
    if n > n_points:
        indices = np.linspace(0, n - 1, n_points, dtype=int)
        pv = pv[indices]
        expected = expected[indices]

    # Convert to -log10 scale
    obs_log = -np.log10(pv)
    exp_log = -np.log10(expected)

    return exp_log, obs_log


def generate_confidence_band(n_total: int, n_points: int = 1000, alpha: float = 0.05) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate 95% confidence band for expected QQ plot.

    Uses order statistics to compute exact confidence intervals.

    Args:
        n_total: Total number of p-values
        n_points: Number of points to return
        alpha: Significance level (0.05 = 95% CI)

    Returns:
        (expected, ci_low, ci_high) in -log10 scale
    """
    expected = np.linspace(1/(n_total+1), n_total/(n_total+1), n_points)

    # Use beta distribution for order statistics
    ci_low = []
    ci_high = []

    for i, e in enumerate(expected):
        k = int(e * n_total)  # approximate rank
        k = max(1, min(k, n_total - 1))

        # Beta distribution for k-th order statistic
        low = stats.beta.ppf(alpha/2, k, n_total - k + 1)
        high = stats.beta.ppf(1 - alpha/2, k, n_total - k + 1)
        ci_low.append(low)
        ci_high.append(high)

    exp_log = -np.log10(expected)
    ci_low_log = -np.log10(np.array(ci_high))  # Note: flip for log scale
    ci_high_log = -np.log10(np.array(ci_low))

    return exp_log, ci_low_log, ci_high_log


def generate_all_qq_data(experiments: List[str] = None, max_samples: int = None) -> Dict[str, pd.DataFrame]:
    """
    Generate QQ plot data for all specified experiments.

    Args:
        experiments: List of experiment IDs (default: all in METHOD_METADATA)
        max_samples: Max samples per experiment (for testing)

    Returns:
        Dictionary with:
        - 'qq_points': DataFrame with expected/observed for each method
        - 'summary': DataFrame with lambda, n_pvalues, metadata per method
        - 'confidence': DataFrame with 95% CI band
    """
    if experiments is None:
        experiments = list(METHOD_METADATA.keys())

    print("=" * 60)
    print("ENHANCED QQ DATA GENERATION")
    print(f"Experiments: {len(experiments)}")
    print("=" * 60)

    qq_data = []
    summary_data = []

    for exp_id in experiments:
        print(f"\nProcessing {exp_id}...")

        # All experiment RESULTS are in stage1_global_dispersion/results/
        # (The stage1_cov_dispersion_bins only contains dispersion parameter files)
        results_dir = GLOBAL_DISP_DIR

        pvalues = load_experiment_pvalues(exp_id, results_dir, max_samples)

        if len(pvalues) == 0:
            print(f"  Skipping {exp_id}: no p-values loaded")
            continue

        print(f"  Loaded {len(pvalues):,} p-values")

        # Compute lambda
        lambda_gc = compute_lambda_gc(pvalues)
        print(f"  Lambda (GC): {lambda_gc:.3f}")

        # Generate QQ points
        exp_log, obs_log = generate_qq_points(pvalues)

        # Store QQ data
        for e, o in zip(exp_log, obs_log):
            qq_data.append({
                'experiment': exp_id,
                'expected_log10': e,
                'observed_log10': o
            })

        # Store summary
        metadata = METHOD_METADATA.get(exp_id, {})
        summary_data.append({
            'experiment': exp_id,
            'name': metadata.get('name', exp_id),
            'dispersion': metadata.get('dispersion', 'Unknown'),
            'fdr_method': metadata.get('fdr', 'Unknown'),
            'estimator': metadata.get('estimator', 'Unknown'),
            'n_pvalues': len(pvalues),
            'lambda_gc': lambda_gc
        })

    # Generate confidence band (using median n_pvalues)
    if summary_data:
        median_n = int(np.median([s['n_pvalues'] for s in summary_data]))
        exp_ci, ci_low, ci_high = generate_confidence_band(median_n)
        confidence_df = pd.DataFrame({
            'expected_log10': exp_ci,
            'ci_low_log10': ci_low,
            'ci_high_log10': ci_high
        })
    else:
        confidence_df = pd.DataFrame()

    return {
        'qq_points': pd.DataFrame(qq_data),
        'summary': pd.DataFrame(summary_data),
        'confidence': confidence_df
    }


def save_qq_data(data: Dict[str, pd.DataFrame], output_dir: Path = None):
    """Save QQ data to TSV files."""
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "data" / "qq_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, df in data.items():
        if not df.empty:
            output_path = output_dir / f"qq_{name}.tsv"
            df.to_csv(output_path, sep='\t', index=False)
            print(f"Saved: {output_path}")


def print_summary_table(summary_df: pd.DataFrame):
    """Print formatted summary table."""
    if summary_df.empty:
        return

    print("\n" + "=" * 80)
    print("METHOD COMPARISON SUMMARY")
    print("=" * 80)
    print(f"{'Experiment':<12} {'Name':<22} {'Dispersion':<18} {'FDR':<10} {'Lambda':>8} {'N P-values':>12}")
    print("-" * 80)

    # Sort by lambda (closest to 1.0 is best)
    summary_df = summary_df.copy()
    summary_df['lambda_distance'] = np.abs(summary_df['lambda_gc'] - 1.0)
    summary_df = summary_df.sort_values('lambda_distance')

    for _, row in summary_df.iterrows():
        # Mark best method
        marker = "*" if row['lambda_distance'] == summary_df['lambda_distance'].min() else " "
        print(f"{row['experiment']:<12} {row['name']:<22} {row['dispersion']:<18} "
              f"{row['fdr_method']:<10} {row['lambda_gc']:>8.3f} {row['n_pvalues']:>12,}{marker}")

    print("-" * 80)
    print("* Best method (lambda closest to 1.0)")
    print()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate enhanced QQ data for method comparison')
    parser.add_argument('--experiments', nargs='+', default=None,
                        help='Specific experiments to process (default: all)')
    parser.add_argument('--max-samples', type=int, default=None,
                        help='Max samples per experiment (for testing)')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory for TSV files')

    args = parser.parse_args()

    # Generate data
    data = generate_all_qq_data(args.experiments, args.max_samples)

    # Print summary
    print_summary_table(data['summary'])

    # Save data
    save_qq_data(data, args.output_dir)

    print("\nDone!")
