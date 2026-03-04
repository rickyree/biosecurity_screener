#!/usr/bin/env python3
"""
TM-align comparison of low-ranking _A candidates vs A_B_b2_noPTM_A.
Focuses on the interface region (residues in chain A within 5 Å of chain B).
"""

import copy, warnings
import numpy as np
from scipy.spatial import KDTree
from Bio.PDB import PDBParser
from tmtools import tm_align

warnings.filterwarnings('ignore')

REFERENCE_COMPLEX = 'binding/A_B_b2_noPTM.pdb'
REFERENCE_A       = 'binding/A_B_b2_noPTM_A.pdb'

CANDIDATES = [
    'binding/2x9m_A.pdb',
    'binding/6cmg_A.pdb',
    'binding/2hle_A.pdb',
    'binding/8wjb_A.pdb',
    'binding/6p72_A.pdb',
    'binding/1v3d_A.pdb',
    'binding/7syy_A.pdb',
    'binding/8vf1_A.pdb',
    'binding/5nop_A.pdb',
    'binding/3gxu_A.pdb',
    'binding/2bba_A.pdb',
    'binding/6vy4_A.pdb',
    'binding/1trn_A.pdb',
    'binding/2vwd_A.pdb',
    'binding/6cmi_A.pdb',
    'binding/1iko_A.pdb',
    'binding/5kv9_A.pdb',
]

RANKS = {
    '2x9m_A': 24, '6cmg_A': 29, '2hle_A': 30, '8wjb_A': 32,
    '6p72_A': 33, '1v3d_A': 37, '7syy_A': 40, '8vf1_A': 42,
    '5nop_A': 44, '3gxu_A': 46, '2bba_A': 55, '6vy4_A': 60,
    '1trn_A': 62, '2vwd_A': 63, '6cmi_A': 65, '1iko_A': 66,
    '5kv9_A': 67,
}

parser = PDBParser(QUIET=True)


def get_ca(structure):
    """Extract Cα coordinates (N,3) and one-letter residue labels."""
    coords, seq = [], []
    for model in structure:
        for chain in model:
            for res in chain:
                if 'CA' in res:
                    coords.append(res['CA'].coord.copy())
                    seq.append(res.get_resname().strip()[:1])  # just a label
    return np.array(coords, dtype=np.float64), ''.join(seq)


def get_iface_resids(complex_path):
    """Return set of (chain_id, res_id) for chain-A residues ≤5 Å from chain B.
    Uses only chain + residue-id tuple (strips structure/model prefix) so it
    matches across separately-loaded PDB files.
    """
    struct = parser.get_structure('cplx', complex_path)
    A_res_ca = {}
    B_coords = []
    for model in struct:
        for chain in model:
            for res in chain:
                if chain.id == 'A' and 'CA' in res:
                    key = (chain.id, res.get_id())   # ('A', (' ', 188, ' '))
                    A_res_ca[key] = res['CA'].coord.copy()
                elif chain.id == 'B':
                    for atom in res:
                        B_coords.append(atom.coord.copy())
    if not B_coords:
        return set()
    tree = KDTree(np.array(B_coords))
    iface = set()
    for key, ca in A_res_ca.items():
        if tree.query_ball_point(ca, r=5.0):
            iface.add(key)
    return iface


def get_ca_ordered(structure):
    """Return per-residue Cα coords and (chain_id, res_id) keys, in chain order."""
    coords, keys = [], []
    for model in structure:
        for chain in model:
            for res in chain:
                if 'CA' in res:
                    coords.append(res['CA'].coord.copy())
                    keys.append((chain.id, res.get_id()))
    return np.array(coords, dtype=np.float64), keys


# ── Build reference ──────────────────────────────────────────────────────────

print("Loading reference A_B_b2_noPTM_A …")
ref_struct = parser.get_structure('ref', REFERENCE_A)
ref_ca, ref_seq = get_ca(ref_struct)
ref_ca_ord, ref_fids = get_ca_ordered(ref_struct)

print("Identifying interface residues from complex …")
iface_fids = get_iface_resids(REFERENCE_COMPLEX)
print(f"  Interface residues on A : {len(iface_fids)}")

# Indices of interface residues in the reference Cα array
iface_idx_ref = np.array([i for i, key in enumerate(ref_fids) if key in iface_fids], dtype=int)
print(f"  Mapped to {len(iface_idx_ref)} Cα positions in reference")

# ── Run TM-align for each candidate ─────────────────────────────────────────

print()
print("=" * 80)
print(f" TM-align vs A_B_b2_noPTM_A  |  interface: {len(iface_idx_ref)} residues ≤5 Å from B")
print("=" * 80)
header = f"{'Rank':>4}  {'Candidate':<14}  {'#res':>5}  {'TM-score':>8}  {'RMSD_all':>8}  {'RMSD_iface':>10}  {'Iface_aln':>10}"
print(header)
print("-" * len(header))

rows = []

for path in CANDIDATES:
    name = path.split('/')[-1].replace('.pdb', '')
    rank = RANKS.get(name, '?')

    cand_struct = parser.get_structure(name, path)
    cand_ca, cand_seq = get_ca(cand_struct)

    if len(cand_ca) < 5:
        print(f"{rank:>4}  {name:<14}  too few residues — skipped")
        continue

    # TM-align: aligns cand onto ref (cand = coords1, ref = coords2)
    result = tm_align(cand_ca, ref_ca, cand_seq, ref_seq)
    tm_score = result.tm_norm_chain2   # normalised by reference length
    rmsd_all = result.rmsd

    # Apply the TM-align rotation to the candidate Cα coords
    R   = result.u          # (3,3) rotation matrix
    t   = result.t          # (3,) translation vector
    cand_aligned = (R @ cand_ca.T).T + t

    # Interface RMSD: find reference interface Cα, then NN in aligned candidate
    ref_iface_ca = ref_ca_ord[iface_idx_ref]
    tree_cand = KDTree(cand_aligned)
    nn_dists, nn_idx = tree_cand.query(ref_iface_ca)

    # Only count pairs within 8 Å (avoid matching to far-away loops)
    cutoff = 8.0
    paired = nn_dists < cutoff
    n_paired = paired.sum()
    if n_paired >= 3:
        iface_rmsd = float(np.sqrt(np.mean(nn_dists[paired]**2)))
    else:
        iface_rmsd = float('inf')

    n_cand = len(cand_ca)
    iface_str = f"{iface_rmsd:.2f}" if np.isfinite(iface_rmsd) else "  —"
    print(f"{rank:>4}  {name:<14}  {n_cand:>5}  {tm_score:>8.4f}  {rmsd_all:>8.2f}  {iface_str:>10}  {n_paired:>5}/{len(iface_idx_ref):<5}")

    rows.append(dict(rank=rank, name=name, n_res=n_cand,
                     tm_score=tm_score, rmsd_all=rmsd_all,
                     iface_rmsd=iface_rmsd, iface_paired=n_paired,
                     iface_total=len(iface_idx_ref)))

# ── Save TSV ─────────────────────────────────────────────────────────────────
tsv_out = 'results/dmasif/tmalign_iface.tsv'
with open(tsv_out, 'w') as f:
    f.write('rank\tcandidate\tn_res\ttm_score\trmsd_all_A\tiface_rmsd_A\tiface_paired\tiface_total\n')
    for r in rows:
        iface_str = f"{r['iface_rmsd']:.3f}" if np.isfinite(r['iface_rmsd']) else 'inf'
        f.write(f"{r['rank']}\t{r['name']}\t{r['n_res']}\t{r['tm_score']:.4f}\t"
                f"{r['rmsd_all']:.3f}\t{iface_str}\t{r['iface_paired']}\t{r['iface_total']}\n")

print(f"\nSaved → {tsv_out}")
print()
print("Legend:")
print("  TM-score   >0.5 = same fold  /  >0.8 = very similar")
print("  RMSD_all   global Cα RMSD after TM-align superposition")
print("  RMSD_iface NN-based Cα RMSD restricted to reference interface residues (<8 Å cutoff)")
print("  Iface_aln  how many interface residues found a near match in candidate")
