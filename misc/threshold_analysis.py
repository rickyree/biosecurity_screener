#!/usr/bin/env python3
"""
First-cluster threshold classifier for dMaSIF scores.

Strategy: scan upward from the lowest score. Find the first position
where the gap before a run of scores is much larger than the gaps
within that run (ratio method). The threshold sits just below the
start of that first dense cluster, so the cluster — and everything
above it — is labelled non-binder.

Usage:
  python threshold_analysis.py --tsv rankings_combined.tsv --score-col combined_score
  python threshold_analysis.py --tsv rankings.tsv --score-col score
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm as sp_norm
from sklearn.mixture import GaussianMixture


# ── Core: find the first dense cluster scanning from the bottom ──────────────

def find_first_cluster_edge(scores, fwd_window=5, ratio_thresh=3.0):
    """Scan upward through sorted scores and find the first gap where:
        gap[i] / mean(gap[i+1 : i+1+fwd_window]) >= ratio_thresh

    That gap marks the transition from scattered points to the first
    dense cluster. The threshold is placed at the midpoint of that gap
    (i.e., just below where the cluster starts).

    Returns
    -------
    threshold : float   — cut-off score
    cluster_start : float — first score inside the first cluster
    trigger_idx : int   — index i where the ratio fired
    ratio : float       — ratio value at that index
    """
    s = np.sort(scores)
    gaps = np.diff(s)
    n = len(gaps)

    for i in range(n - fwd_window):
        fwd_mean = gaps[i + 1: i + 1 + fwd_window].mean()
        if fwd_mean < 1e-10:
            continue
        ratio = gaps[i] / fwd_mean
        if ratio >= ratio_thresh:
            threshold    = float((s[i] + s[i + 1]) / 2)
            cluster_start = float(s[i + 1])
            return threshold, cluster_start, i, ratio

    # Fallback: midpoint of largest gap
    gi = int(np.argmax(gaps))
    return float((s[gi] + s[gi + 1]) / 2), float(s[gi + 1]), gi, float(gaps[gi] / (gaps.mean() + 1e-12))


# ── GMM (2 components) for visualisation only ────────────────────────────────

def fit_gmm(scores):
    gm = GaussianMixture(n_components=2, random_state=42).fit(scores.reshape(-1, 1))
    means = gm.means_.flatten()
    stds  = np.sqrt(gm.covariances_.flatten())
    wts   = gm.weights_
    lo, hi = np.argsort(means)
    return {
        'binder':    {'mu': means[lo], 'sigma': stds[lo], 'w': wts[lo]},
        'nonbinder': {'mu': means[hi], 'sigma': stds[hi], 'w': wts[hi]},
    }, gm


# ── Classify ─────────────────────────────────────────────────────────────────

def classify(df, score_col, threshold):
    df = df.copy()
    df['label'] = np.where(df[score_col] < threshold, 'binder', 'non-binder')
    return df


# ── Plot ─────────────────────────────────────────────────────────────────────

def plot_results(df, score_col, threshold, cluster_start, trigger_idx, components, out_path):
    scores = df[score_col].values
    s_sorted = np.sort(scores)
    xs = np.linspace(s_sorted.min() - 0.1, s_sorted.max() + 0.1, 600)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # ── Left: score distribution + GMM overlay ───────────────────────────────
    ax = axes[0]
    ax.hist(scores, bins=28, density=True, color='#b0c4de', alpha=0.55,
            edgecolor='white', label='scores')

    bn = components['binder']
    nb = components['nonbinder']
    ax.plot(xs, bn['w'] * sp_norm.pdf(xs, bn['mu'], bn['sigma']),
            color='steelblue', lw=1.8, ls='--', alpha=0.8,
            label=f"GMM low  μ={bn['mu']:+.3f}")
    ax.plot(xs, nb['w'] * sp_norm.pdf(xs, nb['mu'], nb['sigma']),
            color='salmon', lw=1.8, ls='--', alpha=0.8,
            label=f"GMM high μ={nb['mu']:+.3f}")

    ax.axvline(threshold, color='crimson', lw=2.5,
               label=f'Threshold  {threshold:+.3f}')
    ax.axvspan(scores.min() - 0.1, threshold,
               color='steelblue', alpha=0.10)
    ax.axvspan(threshold, scores.max() + 0.1,
               color='salmon', alpha=0.10)

    ax.set_xlabel('Score (lower = better)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('Score distribution — first-cluster threshold', fontsize=12)
    ax.legend(fontsize=8.5, loc='upper left')

    # ── Right: waterfall with classification ─────────────────────────────────
    ax2 = axes[1]
    df_sorted = df.sort_values(score_col).reset_index(drop=True)
    colors = ['steelblue' if l == 'binder' else 'salmon'
              for l in df_sorted['label']]
    ax2.scatter(range(len(df_sorted)), df_sorted[score_col],
                c=colors, s=42, zorder=3, alpha=0.85)
    ax2.plot(range(len(df_sorted)), df_sorted[score_col],
             color='gray', lw=0.8, alpha=0.45, zorder=2)

    ax2.axhline(threshold, color='crimson', lw=2.5, ls='--',
                label=f'Threshold  {threshold:+.3f}')
    ax2.axhline(cluster_start, color='darkorange', lw=1.2, ls=':',
                label=f'Cluster start  {cluster_start:+.3f}')

    n_b = (df_sorted['label'] == 'binder').sum()
    n_n = (df_sorted['label'] == 'non-binder').sum()
    mid_b = n_b / 2
    mid_n = n_b + n_n / 2
    ax2.text(mid_b, threshold - 0.07, f'{n_b} binders',
             ha='center', fontsize=9, color='steelblue', fontweight='bold')
    ax2.text(mid_n, threshold + 0.07, f'{n_n} non-binders',
             ha='center', fontsize=9, color='salmon', fontweight='bold')

    # Annotate the triggering gap
    ax2.annotate(f'ratio={components.get("_ratio", 0):.1f}',
                 xy=(trigger_idx + 0.5, (threshold + cluster_start) / 2),
                 fontsize=8, color='crimson', ha='left')

    ax2.set_xlabel('Rank (sorted by score)', fontsize=11)
    ax2.set_ylabel('Score', fontsize=11)
    ax2.set_title('Waterfall — first-cluster classification', fontsize=12)
    ax2.legend(fontsize=8.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"  Figure saved → {out_path}")


# ── Main ─────────────────────────────────────────────────────────────────────

SCORE_COL_CANDIDATES = ['score', 'combined_score']


def detect_score_col(df, explicit=None):
    if explicit:
        if explicit not in df.columns:
            raise SystemExit(f"ERROR: column '{explicit}' not found. Available: {list(df.columns)}")
        return explicit
    for col in SCORE_COL_CANDIDATES:
        if col in df.columns:
            return col
    raise SystemExit(f"ERROR: no score column found. Expected one of {SCORE_COL_CANDIDATES}. Got: {list(df.columns)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tsv',
                    help='Rankings TSV file to classify (score column auto-detected)')
    ap.add_argument('--score-col',  default=None,
                    help='Score column name (auto-detected if omitted)')
    ap.add_argument('--fwd-window', type=int,   default=5,
                    help='Number of forward gaps used to estimate cluster density (default: 5)')
    ap.add_argument('--ratio',      type=float, default=3.0,
                    help='Gap/forward-mean ratio that triggers cluster detection (default: 3.0)')
    ap.add_argument('--out-fig',    default=None)
    ap.add_argument('--out-tsv',    default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep='\t')
    score_col = detect_score_col(df, args.score_col)
    args.score_col = score_col
    df[score_col] = pd.to_numeric(df[score_col], errors='coerce')
    df = df.dropna(subset=[score_col])
    scores = df[score_col].values

    print(f"\n{'='*62}")
    print(f" First-cluster threshold classifier")
    print(f" TSV : {args.tsv}   col: {score_col} (auto-detected)")
    print(f"{'='*62}")
    print(f"  N = {len(scores)}  |  range [{scores.min():+.4f}, {scores.max():+.4f}]")
    print(f"  Parameters: fwd_window={args.fwd_window}, ratio_thresh={args.ratio}")

    threshold, cluster_start, trig_i, ratio = find_first_cluster_edge(
        scores, fwd_window=args.fwd_window, ratio_thresh=args.ratio)

    components, _ = fit_gmm(scores)
    components['_ratio'] = ratio   # pass through for annotation

    print(f"\n  First cluster detected:")
    print(f"    Trigger index   : {trig_i}  (score {np.sort(scores)[trig_i]:+.4f})")
    print(f"    Cluster starts  : {cluster_start:+.4f}")
    print(f"    Gap/fwd ratio   : {ratio:.2f}  (threshold={args.ratio})")
    print(f"    → Threshold set : {threshold:+.4f}  (midpoint of triggering gap)")

    df_out = classify(df, args.score_col, threshold)
    n_b = (df_out['label'] == 'binder').sum()
    n_n = (df_out['label'] == 'non-binder').sum()

    print(f"\n  Classification at {threshold:+.4f}:")
    print(f"    BINDER     : {n_b}")
    print(f"    non-binder : {n_n}")

    print(f"\n  Binders ({n_b}):")
    binders = df_out[df_out['label'] == 'binder'].sort_values(args.score_col)
    for _, row in binders.iterrows():
        print(f"    {row['candidate']:<40} {row[args.score_col]:+.4f}")

    stem    = args.tsv.replace('.tsv', '')
    out_fig = args.out_fig or f"{stem}_first_cluster_threshold.png"
    out_tsv = args.out_tsv or f"{stem}_first_cluster_classified.tsv"

    plot_results(df_out, args.score_col, threshold, cluster_start,
                 trig_i, components, out_fig)
    df_out.to_csv(out_tsv, sep='\t', index=False)
    print(f"  Classified TSV → {out_tsv}")
    print()


if __name__ == '__main__':
    main()
