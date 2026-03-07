import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv('results/to_run_results3.tsv', sep='\t')
df['score'] = df['chem_sim']

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.facecolor': '#ffffff',
    'figure.facecolor': '#ffffff',
    'axes.grid': True,
    'grid.color': '#e2e8f0',
    'grid.linewidth': 0.8,
})

fig, ax = plt.subplots(figsize=(10, 5.5))

sns.histplot(
    data=df, x='score',
    bins=30, edgecolor='white', linewidth=0.5,
    ax=ax, color='#60a5fa'
)

ax.set_xlim(0.7, 0.9)
ax.set_title('ESM2 Chemical Similarity Scores - Cholix Full-length Variants',
             fontsize=13, fontweight='bold', color='#0f172a', pad=12)
ax.set_xlabel('Chemical Similarity Score', fontsize=11, color='#475569')
ax.set_ylabel('Number of sequences', fontsize=11, color='#475569')
ax.tick_params(colors='#475569')

plt.tight_layout()
plt.savefig('cholix_esm2_histogram.png', dpi=180, bbox_inches='tight')
print("Saved cholix_esm2_histogram.png")
