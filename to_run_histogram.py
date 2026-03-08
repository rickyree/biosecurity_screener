import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import seaborn as sns

# ── Data ──────────────────────────────────────────────────────────────────────


filename = 'output.tsv'


df = pd.read_csv(filename, sep='\t')

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



ax.set_xlim(0, 1)
ax.set_title('Aegis scores of Cholix toxin synthetic homologues',
             fontsize=13, fontweight='bold', color='#0f172a', pad=12)
ax.set_xlabel('Normalised Aegis score', fontsize=11, color='#475569')
ax.set_ylabel('Number of sequences', fontsize=11, color='#475569')
ax.tick_params(colors='#475569')

patch_flagged = mpatches.Patch(color='#f87171', label='Flagged by Commec', ec='white')
patch_passed  = mpatches.Patch(color='#60a5fa', label='Passed by Commec',  ec='white')

ax.legend(handles=[patch_passed, patch_flagged],
          frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig('aegis__histogram_.png', dpi=180, bbox_inches='tight')
print("Saved aegis_histogram.png")
