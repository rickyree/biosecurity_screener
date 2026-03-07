import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import seaborn as sns

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv('to_run_results_all.tsv', sep='\t')

with open('to_run_results_bsspassed.tsv') as f:
    passed = set(line.strip() for line in f if line.strip())

df['Commec'] = df['candidate'].apply(
    lambda x: 'Passed by Commec' if x in passed else 'Flagged by Commec'
)

mn, mx = df['chem_sim'].min(), df['chem_sim'].max()
df['score'] = (df['chem_sim'] - mn) / (mx - mn)
df['AegisLabel'] = df['score'].apply(lambda x: 'Flagged' if x >= 0.3 else 'Passed')

THRESHOLD = 0.3

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

fig = plt.figure(figsize=(14, 5))
gs = gridspec.GridSpec(1, 2, width_ratios=[2.2, 1], wspace=0.35)

# ── Left: Histogram ───────────────────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0])

palette = {'Flagged by Commec': '#f87171', 'Passed by Commec': '#60a5fa'}
sns.histplot(
    data=df, x='score', hue='Commec',
    multiple='stack', palette=palette,
    bins=30, edgecolor='white', linewidth=0.5,
    ax=ax1,
)

ax1.axvline(THRESHOLD, color='#0f172a', linestyle='--', linewidth=1.8)
ax1.text(THRESHOLD + 0.01, ax1.get_ylim()[1] * 0.95,
         f'Aegis threshold ({THRESHOLD})',
         color='#0f172a', fontsize=9, va='top')

ax1.set_xlim(0, 1)
ax1.set_title('Aegis scores of Cholix toxin synthetic homologues',
              fontsize=12, fontweight='bold', color='#0f172a', pad=10)
ax1.set_xlabel('Normalised Aegis score', fontsize=10, color='#475569')
ax1.set_ylabel('Number of sequences', fontsize=10, color='#475569')
ax1.tick_params(colors='#475569')

patch_flagged = mpatches.Patch(color='#f87171', label='Flagged by Commec', ec='white')
patch_passed  = mpatches.Patch(color='#60a5fa', label='Passed by Commec',  ec='white')
line_thresh   = mlines.Line2D([], [], color='#0f172a', linestyle='--',
                               linewidth=1.8, label=f'Aegis threshold ({THRESHOLD})')
ax1.legend(handles=[line_thresh, patch_passed, patch_flagged],
           frameon=False, fontsize=9)

# ── Right: Confusion matrix ───────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[1])

labels = ['Flagged', 'Passed']
cm = pd.crosstab(df['AegisLabel'], df['Commec'].apply(
    lambda x: 'Flagged' if x == 'Flagged by Commec' else 'Passed'
)).reindex(index=labels, columns=labels, fill_value=0)

cmap = mcolors.LinearSegmentedColormap.from_list('aegis', ['#f0f9ff', '#0284c7'])
sns.heatmap(
    cm, annot=True, fmt='d', cmap=cmap, cbar=False,
    linewidths=2, linecolor='white',
    annot_kws={'size': 22, 'weight': 'bold', 'color': '#0f172a'},
    ax=ax2,
)

ax2.set_title('')
ax2.set_xlabel('Commec', fontsize=14, color='#475569', labelpad=10)
ax2.set_ylabel('Aegis', fontsize=14, color='#475569', labelpad=10)
ax2.tick_params(colors='#475569', length=0)
ax2.set_xticklabels([u'\u2717', u'\u2713'], fontsize=20)
ax2.set_yticklabels([u'\u2717', u'\u2713'], fontsize=20, rotation=0)

for spine in ax2.spines.values():
    spine.set_visible(False)

# ── Save ──────────────────────────────────────────────────────────────────────
plt.savefig('aegis_cholix_results.png', dpi=180, bbox_inches='tight')
print("Saved aegis_cholix_results.png")
