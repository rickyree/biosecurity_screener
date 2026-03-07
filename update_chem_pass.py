#!/usr/bin/env python3
"""
Update the chem_pass column in to_run_results_all.tsv based on new threshold:
chem_pass = True if chem_sim < 0.755, False otherwise
"""
import pandas as pd

# Read the data
df = pd.read_csv('results/to_run_results_all.tsv', sep='\t')

print(f"Original data: {len(df)} rows")
print(f"Original chem_pass counts:")
print(df['chem_pass'].value_counts())

# Update chem_pass column with new logic: True if chem_sim < 0.755
NEW_THRESHOLD = 0.755
df['chem_pass'] = df['chem_sim'] < NEW_THRESHOLD

print(f"\nAfter update (threshold < {NEW_THRESHOLD}):")
print(df['chem_pass'].value_counts())

# Update pass_filter column to reflect new chem_pass logic
df['pass_filter'] = df['contact_pass'] & df['chem_pass']

print(f"\nUpdated pass_filter counts:")
print(df['pass_filter'].value_counts())

# Save back to the same file
df.to_csv('results/to_run_results_all.tsv', sep='\t', index=False)

print(f"\nUpdated file saved: results/to_run_results_all.tsv")
print(f"New logic: chem_pass = True if chem_sim < {NEW_THRESHOLD}")