"""
Summarise USalign results for cholix_fullprotein_t10 mutconserved sequences.

Reads:  results/usalign/mutconserved_t10/metrics.txt
Writes: results/usalign/mutconserved_t10/summary_per_comparison.csv
        results/usalign/mutconserved_t10/summary_per_mobile.csv
Prints: table to stdout
"""

import re, csv
from collections import defaultdict

METRICS       = "results/usalign/mutconserved_t10/metrics.txt"
OUT_DIR       = "results/usalign/mutconserved_t10"
THRESH_FLAG   = 0.7
THRESH_BORDER = 0.5

# ---------- parse metrics.txt ----------
text   = open(METRICS).read()
blocks = re.split(r"=== (.+?) ===", text)

rows = []
for i in range(1, len(blocks), 2):
    job_name = blocks[i].strip()
    content  = blocks[i + 1]

    def get(pattern):
        m = re.search(pattern, content)
        return float(m.group(1)) if m else None

    tm1   = get(r"tm_score_1: ([0-9.]+)")
    tm2   = get(r"tm_score_2: ([0-9.]+)")
    rmsd  = get(r"rmsd: ([0-9.]+)")
    seqid = get(r"sequence_identity: ([0-9.]+)")
    alen  = get(r"aligned_length: ([0-9]+)")

    if tm1 is None or tm2 is None:
        print(f"  WARNING: could not parse {job_name}")
        continue

    m = re.match(r"(cholix_fullprotein_t10_mutconserved_\d+) vs (.+)", job_name)
    mobile = m.group(1) if m else job_name
    target = m.group(2) if m else ""

    rows.append({
        "mobile":        mobile,
        "target":        target,
        "tm_score_1":    tm1,
        "tm_score_2":    tm2,
        "rmsd":          rmsd,
        "seq_id":        seqid,
        "aligned_len":   alen,
        "classification": "flagged"    if tm2 >= THRESH_FLAG   else
                          "borderline" if tm2 >= THRESH_BORDER else "pass",
    })

print(f"Parsed {len(rows)} comparisons\n")

# ---------- write per-comparison CSV ----------
comp_csv = f"{OUT_DIR}/summary_per_comparison.csv"
with open(comp_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"Written: {comp_csv}")

# ---------- per-mobile summary ----------
mobile_stats = defaultdict(lambda: {
    "mobile": "", "max_tm2": 0.0, "best_target": "",
    "n_flagged": 0, "n_borderline": 0, "n_targets": 0,
})

for r in rows:
    s = mobile_stats[r["mobile"]]
    s["mobile"]    = r["mobile"]
    s["n_targets"] += 1
    if r["classification"] == "flagged":
        s["n_flagged"] += 1
    elif r["classification"] == "borderline":
        s["n_borderline"] += 1
    if r["tm_score_2"] > s["max_tm2"]:
        s["max_tm2"]     = r["tm_score_2"]
        s["best_target"] = r["target"]

mob_rows = sorted(mobile_stats.values(), key=lambda x: -x["max_tm2"])

mob_csv = f"{OUT_DIR}/summary_per_mobile.csv"
with open(mob_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["mobile", "max_tm2", "best_target", "n_flagged", "n_borderline", "n_targets"])
    w.writeheader()
    w.writerows(mob_rows)
print(f"Written: {mob_csv}\n")

# ---------- printed summary ----------
def classify(tm2):
    if tm2 >= THRESH_FLAG:   return "FLAGGED"
    if tm2 >= THRESH_BORDER: return "BORDERLINE"
    return ""

flagged    = [e for e in mob_rows if e["max_tm2"] >= THRESH_FLAG]
borderline = [e for e in mob_rows if THRESH_BORDER <= e["max_tm2"] < THRESH_FLAG]

print(f"{'='*75}")
print(f"  cholix_fullprotein_t10 mutconserved  —  "
      f"{len(flagged)}/{len(mob_rows)} flagged (>={THRESH_FLAG}), "
      f"{len(borderline)} borderline ({THRESH_BORDER}–{THRESH_FLAG})")
print(f"{'='*75}")
print(f"  {'mobile':<45} {'max_tm2':>8}  {'hits/3':>6}  best_target")
print(f"  {'-'*44} {'-'*8}  {'-'*6}  {'-'*30}")
for e in mob_rows:
    label = classify(e["max_tm2"])
    suffix = f" *** {label}" if label == "FLAGGED" else f" *  {label}" if label else ""
    print(f"  {e['mobile']:<45} {e['max_tm2']:>8.4f}  {e['n_flagged']:>3}/{e['n_targets']:<2}  {e['best_target']}{suffix}")

print(f"{'='*75}")
print(f"  Flagged    : {len(flagged)}/{len(mob_rows)}")
print(f"  Borderline : {len(borderline)}/{len(mob_rows)}")
print(f"  Pass       : {len(mob_rows) - len(flagged) - len(borderline)}/{len(mob_rows)}")
print(f"{'='*75}")
