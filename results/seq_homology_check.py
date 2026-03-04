#!/usr/bin/env python3
"""
Pairwise sequence alignment of ALL binding/ structures vs A_B_b2_noPTM_A.
Purpose: check whether dMaSIF rankings correlate with sequence homology,
which would compromise the sequence-independence claim.
"""

import glob, warnings
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner

warnings.filterwarnings('ignore')

REFERENCE = 'binding/A_B_b2_noPTM_A.pdb'

# Load ranks from rankings TSV (name -> rank); unranked files get rank=None
rank_map = {}
try:
    with open('results/dmasif/rankings.tsv') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                rank_map[parts[1]] = int(parts[0])
except FileNotFoundError:
    pass

# All PDB files in binding/, sorted; skip the reference itself
all_pdbs = sorted(glob.glob('binding/*.pdb'))
CANDIDATES = [
    (p.replace('binding/', '').replace('.pdb', ''),
     rank_map.get(p.replace('binding/', '').replace('.pdb', ''), None))
    for p in all_pdbs
]

parser  = PDBParser(QUIET=True)
builder = PPBuilder()

aligner = PairwiseAligner()
aligner.mode             = 'global'
aligner.substitution_matrix = __import__('Bio.Align.substitution_matrices',
                                          fromlist=['load']).load('BLOSUM62')
aligner.open_gap_score   = -10
aligner.extend_gap_score = -0.5


def extract_seq(path, use_chain=None):
    """Return the longest polypeptide sequence from the structure (chain A preferred)."""
    struct = parser.get_structure('s', path)
    seqs = {}
    for model in struct:
        for chain in model:
            pp_list = builder.build_peptides(chain)
            if not pp_list:
                continue
            full = ''.join(str(pp.get_sequence()) for pp in pp_list)
            seqs[chain.id] = full
    if not seqs:
        return '', '?'
    if use_chain and use_chain in seqs:
        cid = use_chain
    elif 'A' in seqs:
        cid = 'A'
    else:
        cid = max(seqs, key=lambda k: len(seqs[k]))
    return seqs[cid], cid


def aln_stats(alignment, ref_seq, query_seq):
    """Return (pct_identity, aln_coverage) from a PairwiseAligner Alignment object.

    Uses alignment.indices — each column gives (ref_pos, query_pos) for aligned
    (non-gap) positions only.  Gaps are simply absent from the index arrays.
    """
    ref_idx, qry_idx = alignment.indices          # shape (2, n_aligned)
    n_aligned = len(ref_idx)
    if n_aligned == 0:
        return 0.0, 0.0
    matched = sum(ref_seq[i] == query_seq[j] for i, j in zip(ref_idx, qry_idx))
    pid = matched / n_aligned * 100
    cov = n_aligned / len(ref_seq) * 100          # fraction of ref covered
    return pid, cov


# ── Extract reference sequence ───────────────────────────────────────────────

ref_seq, ref_chain = extract_seq(REFERENCE)
print(f"Reference: A_B_b2_noPTM_A  chain={ref_chain}  len={len(ref_seq)}")
print(f"  {ref_seq[:60]}{'…' if len(ref_seq)>60 else ''}")
print()

# ── Align every candidate ────────────────────────────────────────────────────

results = []
for name, rank in CANDIDATES:
    path = f'binding/{name}.pdb'
    seq, chain = extract_seq(path)
    if not seq:
        print(f"  rank {rank:>2}  {name:<20}  NO SEQUENCE extracted")
        results.append(dict(rank=rank, name=name, chain='?', seq_len=0,
                            pct_id=float('nan'), aln_score=float('nan'),
                            aln_cov=float('nan'), ref_len=len(ref_seq)))
        continue

    alignments = aligner.align(ref_seq, seq)
    best = alignments[0]
    pid, cov = aln_stats(best, ref_seq, seq)

    results.append(dict(rank=rank, name=name, chain=chain, seq_len=len(seq),
                        pct_id=pid, aln_score=float(best.score),
                        aln_cov=cov, ref_len=len(ref_seq)))

# ── Print summary ────────────────────────────────────────────────────────────

def homology_label(pid):
    if np.isnan(pid):    return 'no_sequence'
    if pid >= 90:        return 'near-identical'
    if pid >= 40:        return 'homolog'
    if pid >= 25:        return 'remote homolog'
    return 'non-homolog'

print("=" * 88)
print(" Pairwise sequence alignment vs A_B_b2_noPTM_A  (BLOSUM62, global)")
print("=" * 88)
print(f"  {'Rank':>4}  {'Candidate':<24}  {'Ch':>2}  {'Len':>5}  "
      f"{'%Identity':>10}  {'Cov%':>6}  {'Score':>10}  Homology")
print(f"  {'-'*4}  {'-'*24}  {'-'*2}  {'-'*5}  "
      f"{'-'*10}  {'-'*6}  {'-'*10}  {'-'*15}")

for r in sorted(results, key=lambda x: (x['rank'] is None, x['rank'] or 9999)):
    rank_s = f"{r['rank']:>4}" if r['rank'] is not None else '   —'
    pid_s  = f"{r['pct_id']:>9.1f}%" if not np.isnan(r['pct_id']) else '        —'
    cov_s  = f"{r['aln_cov']:>5.1f}%" if not np.isnan(r['aln_cov']) else '     —'
    sc_s   = f"{r['aln_score']:>10.1f}" if not np.isnan(r['aln_score']) else '         —'
    print(f"  {rank_s}  {r['name']:<24}  {r['chain']:>2}  {r['seq_len']:>5}  "
          f"{pid_s}  {cov_s}  {sc_s}  {homology_label(r['pct_id'])}")

# ── Correlation with rank (ranked files only) ────────────────────────────────

from scipy.stats import spearmanr
ranked = [r for r in results if r['rank'] is not None and not np.isnan(r['pct_id'])]
pids  = np.array([r['pct_id'] for r in ranked])
ranks = np.array([r['rank']   for r in ranked])
rho, pval = spearmanr(ranks, pids)
print()
print(f"  Spearman ρ (dMaSIF rank vs %identity, n={len(ranked)}): {rho:+.3f}  (p={pval:.3e})")
if abs(rho) > 0.5 and pval < 0.05:
    print("  ⚠  Significant correlation — sequence homology may be inflating rankings.")
else:
    print("  ✓  No significant rank–identity correlation.")

# ── Save TSV ─────────────────────────────────────────────────────────────────

tsv_out = 'results/dmasif/seq_homology_top15.tsv'
with open(tsv_out, 'w') as f:
    f.write('rank\tcandidate\tchain\tseq_len\tref_len\tpct_identity\taln_cov_pct\taln_score\thomology_class\n')
    for r in sorted(results, key=lambda x: (x['rank'] is None, x['rank'] or 9999)):
        rank_s = str(r['rank']) if r['rank'] is not None else ''
        pid_s  = f"{r['pct_id']:.2f}"    if not np.isnan(r['pct_id'])    else ''
        cov_s  = f"{r['aln_cov']:.2f}"   if not np.isnan(r['aln_cov'])   else ''
        sc_s   = f"{r['aln_score']:.1f}" if not np.isnan(r['aln_score']) else ''
        f.write(f"{rank_s}\t{r['name']}\t{r['chain']}\t{r['seq_len']}\t{r['ref_len']}\t"
                f"{pid_s}\t{cov_s}\t{sc_s}\t{homology_label(r['pct_id'])}\n")

print(f"\nSaved → {tsv_out}  ({len(results)} entries)")
