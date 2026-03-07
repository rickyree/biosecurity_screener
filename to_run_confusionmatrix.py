import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv('to_run_results_all.tsv', sep='\t')

with open('to_run_results_bsspassed.tsv') as f:
    passed = set(line.strip() for line in f if line.strip())

df['Commec'] = df['candidate'].apply(
    lambda x: 'Passed' if x in passed else 'Flagged'
)

mn, mx = df['chem_sim'].min(), df['chem_sim'].max()
df['score'] = (df['chem_sim'] - mn) / (mx - mn)
df['Aegis'] = df['score'].apply(lambda x: 'Flagged' if x >= 0.3 else 'Passed')

labels = ['Flagged', 'Passed']
cm = pd.crosstab(df['Aegis'], df['Commec']).reindex(
    index=labels, columns=labels, fill_value=0
)

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'figure.facecolor': '#ffffff',
    'axes.facecolor': '#ffffff',
})

cmap = mcolors.LinearSegmentedColormap.from_list(
    'aegis', ['#f0f9ff', '#0284c7']
)

fig, ax = plt.subplots(figsize=(6, 5))

sns.heatmap(
    cm, annot=True, fmt='d', cmap=cmap, cbar=False,
    linewidths=2, linecolor='white',
    annot_kws={'size': 22, 'weight': 'bold', 'color': '#0f172a'},
    ax=ax,
)

ax.set_title('')
ax.set_xlabel('Commec', fontsize=16, color='#475569', labelpad=12)
ax.set_ylabel('Aegis', fontsize=16, color='#475569', labelpad=12)
ax.tick_params(colors='#475569', length=0)
ax.set_xticklabels(['✗', '✓'], fontsize=22)
ax.set_yticklabels(['✗', '✓'], fontsize=22, rotation=0)

for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig('aegis_commec_confusion_matrix.png', dpi=180, bbox_inches='tight')
print("Saved aegis_commec_confusion_matrix.png")