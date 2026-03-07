import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

df = pd.read_csv('to_run_results_all.tsv', sep='\t')

# Extract group labels
df['mut_type'] = df['candidate'].str.extract(r'(mutoutside|mutconserved)')
df['temp'] = df['candidate'].str.extract(r'(lowtemp|hightemp)')
df['group'] = df['mut_type'] + '\n' + df['temp']

group_order = [
    'mutoutside\nlowtemp',
    'mutoutside\nhightemp',
    'mutconserved\nlowtemp',
    'mutconserved\nhightemp',
]

data = [df[df['group'] == g]['chem_sim'].values for g in group_order]

fig, ax = plt.subplots(figsize=(8, 6))

colors = ['#4C9BE8', '#E8834C', '#4CE8A0', '#E84C9B']
bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                medianprops=dict(color='black', linewidth=2))

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Overlay individual points
for i, (d, color) in enumerate(zip(data, colors), start=1):
    jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(d))
    ax.scatter(i + jitter, d, color=color, alpha=0.6, s=20, zorder=3)

labels = ['mutoutside\nlowtemp', 'mutoutside\nhightemp', 'mutconserved\nlowtemp', 'mutconserved\nhightemp']
ax.set_xticks(range(1, len(labels) + 1))
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel('Chemical Similarity (chem_sim)', fontsize=12)
ax.set_title('Chemical Similarity by Mutation Type and Temperature', fontsize=13)
ax.set_ylim(0.7, 0.9)
ax.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('chem_sim_boxplot.png', dpi=150)
print("Saved chem_sim_boxplot.png")
