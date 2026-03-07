#!/usr/bin/env python3
"""
Create box plots for ESM2 prefiltering results comparing mutoutside vs mutconserved sequences.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the data
df = pd.read_csv('to_run_results2.tsv', sep='\t')

# Extract sequence type from candidate names
df['sequence_type'] = df['candidate'].apply(lambda x: 'Mutoutside' if 'mutoutside' in x else 'Mutconserved')

# Set up the plotting style
plt.style.use('default')
sns.set_palette("Set2")

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Define colors for consistency
colors = {'Mutoutside': '#66c2a5', 'Mutconserved': '#fc8d62'}

# Box plot for Contact Score
bp1 = ax1.boxplot([df[df['sequence_type'] == 'Mutoutside']['contact_score'],
                   df[df['sequence_type'] == 'Mutconserved']['contact_score']],
                  labels=['Mutoutside', 'Mutconserved'],
                  patch_artist=True,
                  boxprops=dict(linewidth=1.5),
                  medianprops=dict(color='black', linewidth=2),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5))

# Color the boxes
bp1['boxes'][0].set_facecolor(colors['Mutoutside'])
bp1['boxes'][1].set_facecolor(colors['Mutconserved'])

# Add threshold line for contact score
ref_contact = 0.0398
thresh_contact = ref_contact * 0.5
ax1.axhline(y=thresh_contact, color='red', linestyle='--', alpha=0.7, linewidth=2,
            label=f'Threshold ({thresh_contact:.4f})')
ax1.axhline(y=ref_contact, color='blue', linestyle=':', alpha=0.7, linewidth=2,
            label=f'Reference ({ref_contact:.4f})')

ax1.set_ylabel('Contact Score', fontsize=12, fontweight='bold')
ax1.set_title('ESM2 Contact Score', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)

# Box plot for Chemistry Similarity
bp2 = ax2.boxplot([df[df['sequence_type'] == 'Mutoutside']['chem_sim'],
                   df[df['sequence_type'] == 'Mutconserved']['chem_sim']],
                  labels=['Mutoutside', 'Mutconserved'],
                  patch_artist=True,
                  boxprops=dict(linewidth=1.5),
                  medianprops=dict(color='black', linewidth=2),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5))

# Color the boxes
bp2['boxes'][0].set_facecolor(colors['Mutoutside'])
bp2['boxes'][1].set_facecolor(colors['Mutconserved'])

# Add threshold line for chemistry similarity
thresh_chem = 0.80
ax2.axhline(y=thresh_chem, color='red', linestyle='--', alpha=0.7, linewidth=2,
            label=f'Threshold ({thresh_chem:.2f})')

ax2.set_ylabel('Chemistry Similarity', fontsize=12, fontweight='bold')
ax2.set_title('ESM2 Chemistry Similarity', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10)

# Add statistics annotations
def add_stats(ax, data1, data2, y_max):
    """Add mean and median statistics as text"""
    stats1 = f"n={len(data1)}\nμ={np.mean(data1):.4f}\nm={np.median(data1):.4f}"
    stats2 = f"n={len(data2)}\nμ={np.mean(data2):.4f}\nm={np.median(data2):.4f}"

    ax.text(0.7, y_max * 0.95, stats1, transform=ax.transData,
            bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.7),
            fontsize=9, ha='center')
    ax.text(1.7, y_max * 0.95, stats2, transform=ax.transData,
            bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.7),
            fontsize=9, ha='center')

# Add statistics for contact scores
contact_out = df[df['sequence_type'] == 'Mutoutside']['contact_score']
contact_con = df[df['sequence_type'] == 'Mutconserved']['contact_score']
add_stats(ax1, contact_out, contact_con, ax1.get_ylim()[1])

# Add statistics for chemistry similarity
chem_out = df[df['sequence_type'] == 'Mutoutside']['chem_sim']
chem_con = df[df['sequence_type'] == 'Mutconserved']['chem_sim']
add_stats(ax2, chem_out, chem_con, ax2.get_ylim()[1])

# Overall formatting
plt.suptitle('ESM2 Prefiltering Results: ProteinMPNN Cholix Sequences',
             fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()

# Add pass/fail statistics
total_out = len(df[df['sequence_type'] == 'Mutoutside'])
pass_out = len(df[(df['sequence_type'] == 'Mutoutside') & (df['pass_filter'] == True)])
total_con = len(df[df['sequence_type'] == 'Mutconserved'])
pass_con = len(df[(df['sequence_type'] == 'Mutconserved') & (df['pass_filter'] == True)])

fig.text(0.5, 0.02, f'Filter Results: Mutoutside {pass_out}/{total_out} passed ({pass_out/total_out*100:.1f}%), '
                   f'Mutconserved {pass_con}/{total_con} passed ({pass_con/total_con*100:.1f}%)',
         ha='center', fontsize=11, style='italic')

# Save the plot
plt.savefig('boxplot_esm2_prefilter_results.png', dpi=300, bbox_inches='tight')
plt.show()

print("Box plot saved as 'boxplot_esm2_prefilter_results.png'")
print(f"\nSummary:")
print(f"Mutoutside sequences: {pass_out}/{total_out} passed ({pass_out/total_out*100:.1f}%)")
print(f"Mutconserved sequences: {pass_con}/{total_con} passed ({pass_con/total_con*100:.1f}%)")
