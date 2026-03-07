import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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

THRESHOLD = 0.3

# ── Style ─────────────────────────────────────────────────────────────────────
palette = {'Flagged by Commec': '#f87171', 'Passed by Commec': '#60a5fa'}

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
    data=df, x='score', hue='Commec',
    multiple='stack', palette=palette,
    bins=30, edgecolor='white', linewidth=0.5,
    ax=ax,
)

ax.axvline(THRESHOLD, color='#0f172a', linestyle='--', linewidth=1.8)
ax.text(THRESHOLD + 0.01, ax.get_ylim()[1] * 0.95,
        f'Aegis threshold ({THRESHOLD})',
        color='#0f172a', fontsize=9, va='top')

ax.set_xlim(0, 1)
ax.set_title('Aegis scores of Cholix toxin synthetic homologues',
             fontsize=13, fontweight='bold', color='#0f172a', pad=12)
ax.set_xlabel('Normalised Aegis score', fontsize=11, color='#475569')
ax.set_ylabel('Number of sequences', fontsize=11, color='#475569')
ax.tick_params(colors='#475569')

patch_flagged = mpatches.Patch(color='#f87171', label='Flagged by Commec', ec='white')
patch_passed  = mpatches.Patch(color='#60a5fa', label='Passed by Commec',  ec='white')
line_thresh   = mlines.Line2D([], [], color='#0f172a', linestyle='--',
                               linewidth=1.8, label=f'Aegis threshold ({THRESHOLD})')
ax.legend(handles=[line_thresh, patch_passed, patch_flagged],
          frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig('aegis_cholix_histogram_updated.png', dpi=180, bbox_inches='tight')
print("Saved aegis_cholix_histogram_updated.png")