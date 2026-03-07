import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, classification_report
import numpy as np

# ── Data loading ──────────────────────────────────────────────────────────────
df = pd.read_csv('results/to_run_results_all.tsv', sep='\t')

# Load Commec pass list
with open('to_run_results_bsspassed.tsv') as f:
    passed = set(line.strip() for line in f if line.strip())

# Create binary classifications
df['commec_pass'] = df['candidate'].apply(lambda x: x in passed)
df['esm2_pass'] = df['chem_pass'].fillna(False)

# ── Confusion Matrix ──────────────────────────────────────────────────────────
y_true = df['commec_pass']  # Commec as ground truth
y_pred = df['esm2_pass']    # ESM2 as prediction

cm = confusion_matrix(y_true, y_pred)
labels = ['Flagged by Commec', 'Passed by Commec']
pred_labels = ['Failed ESM2', 'Passed ESM2']

# ── Visualization ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.facecolor': '#ffffff',
    'figure.facecolor': '#ffffff',
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Count matrix
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=pred_labels, yticklabels=labels,
            cbar_kws={'label': 'Count'}, ax=ax1)
ax1.set_title('Confusion Matrix: ESM2 vs Commec\n(Counts)',
              fontsize=13, fontweight='bold', pad=15)
ax1.set_xlabel('Aegis score', fontsize=11)
ax1.set_ylabel('Commec score', fontsize=11)

# Percentage matrix
cm_pct = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100
sns.heatmap(cm_pct, annot=True, fmt='.1f', cmap='Blues',
            xticklabels=pred_labels, yticklabels=labels,
            cbar_kws={'label': 'Percentage (%)'}, ax=ax2)
ax2.set_title('Confusion Matrix: ESM2 vs Commec\n(Row Percentages)',
              fontsize=13, fontweight='bold', pad=15)
ax2.set_xlabel('Aegis score', fontsize=11)
ax2.set_ylabel('Commec score', fontsize=11)

plt.tight_layout()
plt.savefig('cholix_esm2_commec_confusion_matrix.png', dpi=180, bbox_inches='tight')
print("Saved cholix_esm2_commec_confusion_matrix.png")

# ── Classification Report ─────────────────────────────────────────────────────
print("\n" + "="*60)
print("CLASSIFICATION REPORT (ESM2 predicting Commec)")
print("="*60)
print(classification_report(y_true, y_pred,
                           target_names=['Flagged', 'Passed'],
                           digits=3))

# ── Detailed Breakdown ────────────────────────────────────────────────────────
tn, fp, fn, tp = cm.ravel()
total = len(df)

print(f"\nDETAILED BREAKDOWN ({total} total sequences)")
print("-" * 60)
print(f"True Negatives  (TN): {tn:3d} | Flagged by both ESM2 and Commec")
print(f"False Positives (FP): {fp:3d} | Passed ESM2, flagged by Commec")
print(f"False Negatives (FN): {fn:3d} | Failed ESM2, passed by Commec")
print(f"True Positives  (TP): {tp:3d} | Passed by both ESM2 and Commec")

# ── Agreement Analysis ─────────────────────────────────────────────────────────
agreement = (tp + tn) / total
disagreement = (fp + fn) / total

print(f"\nAGREEMENT ANALYSIS")
print("-" * 60)
print(f"Agreement:     {agreement:.3f} ({agreement*100:.1f}%) | {tp + tn}/{total}")
print(f"Disagreement:  {disagreement:.3f} ({disagreement*100:.1f}%) | {fp + fn}/{total}")

if fp > 0:
    print(f"\nESM2 more permissive: {fp} sequences passed ESM2 but flagged by Commec")
if fn > 0:
    print(f"ESM2 more restrictive: {fn} sequences failed ESM2 but passed by Commec")

# ── Sample Disagreements ──────────────────────────────────────────────────────
if fp > 0:
    fp_samples = df[(df['commec_pass'] == False) & (df['esm2_pass'] == True)]
    print(f"\nSAMPLE: ESM2 passed but Commec flagged ({min(5, len(fp_samples))} of {len(fp_samples)}):")
    for _, row in fp_samples.head(5).iterrows():
        cs = row.get('contact_score', 'N/A')
        chem = row.get('chem_sim', 'N/A')
        print(f"  {row['candidate']:<45} | contact={cs} chem_sim={chem}")

if fn > 0:
    fn_samples = df[(df['commec_pass'] == True) & (df['esm2_pass'] == False)]
    print(f"\nSAMPLE: ESM2 failed but Commec passed ({min(5, len(fn_samples))} of {len(fn_samples)}):")
    for _, row in fn_samples.head(5).iterrows():
        cs = row.get('contact_score', 'N/A')
        chem = row.get('chem_sim', 'N/A')
        print(f"  {row['candidate']:<45} | contact={cs} chem_sim={chem}")
