"""
Plot TM-score distribution for cholix t1.5-passed sequences.
Reads: results/usalign/cholix_t1.5_passed/summary_per_mobile.csv
Saves: results/usalign/cholix_t1.5_passed/tm_score_distributions.png
"""

import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

CSV_PATH = "results/usalign/cholix_t1.5_passed/summary_per_mobile.csv"
OUT_PATH = "results/usalign/cholix_t1.5_passed/tm_score_distributions.png"

THRESH_FLAG   = 0.7
THRESH_BORDER = 0.5

# Load data
mobiles = []
vals    = []
with open(CSV_PATH) as fh:
    for r in csv.DictReader(fh):
        mobiles.append(r["mobile"])
        vals.append(float(r["max_tm2"]))

vals    = np.array(vals)
mobiles = np.array(mobiles)

fig, ax = plt.subplots(figsize=(9, 3.5))

bins = np.linspace(0, 1, 21)
ax.hist(vals, bins=bins, color="#1f77b4", alpha=0.8, edgecolor="white", linewidth=0.6)

ax.axvline(THRESH_FLAG,   color="red",    linestyle="--", linewidth=1.2, alpha=0.7,
           label=f"Flagged threshold ({THRESH_FLAG})")
ax.axvline(THRESH_BORDER, color="orange", linestyle="--", linewidth=1.2, alpha=0.7,
           label=f"Borderline threshold ({THRESH_BORDER})")


n_flagged = int(np.sum(vals >= THRESH_FLAG))
n_border  = int(np.sum((vals >= THRESH_BORDER) & (vals < THRESH_FLAG)))

ax.set_xlabel("Max TM-score", fontsize=11)
ax.set_ylabel(f"Count (n={len(vals)})", fontsize=10)
ax.set_xlim(0, 1)
ax.set_ylim(bottom=0, top=ax.get_ylim()[1] * 1.2)
ax.set_title(
    f"TM-score distribution of Cholix toxic sequences (t1.5) flagged by Aegis and passed by Commec\n",
    fontsize=11, fontweight="bold",
)
ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=9, frameon=False)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT_PATH}")
