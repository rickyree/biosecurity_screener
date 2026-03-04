#!/usr/bin/env python3
"""
Compare surface fingerprint properties at the interface region between:
  - Group TOP  : ranks 1–15 (top-performing binders)
  - Group G1   : ranks >16, same fold as reference (TM-score > 0.8), low dMaSIF score

For each structure:
  1. Compute SASA on chain A only.
  2. TM-align chain-A Cα to reference (A_B_b2_noPTM_A).
  3. Transform surface atom coords into reference frame.
  4. Collect atoms within 8 Å of the reference interface centroid.
  5. Record their 9-dim descriptors.

Then compare per-descriptor distributions between groups with a Mann-Whitney U test.
"""

import warnings
import numpy as np
import copy
from scipy.spatial import KDTree
from scipy.stats import mannwhitneyu
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.SASA import ShrakeRupley
from tmtools import tm_align

warnings.filterwarnings('ignore')
np.random.seed(42)

# ── Descriptor tables (same as dmasif_workflow) ──────────────────────────────

KD     = {'ILE':4.5,'VAL':4.2,'LEU':3.8,'PHE':2.8,'CYS':2.5,'MET':1.9,'ALA':1.8,
          'GLY':-0.4,'THR':-0.7,'SER':-0.8,'TRP':-0.9,'TYR':-1.3,'PRO':-1.6,
          'HIS':-3.2,'GLU':-3.5,'GLN':-3.5,'ASP':-3.5,'ASN':-3.5,'LYS':-3.9,'ARG':-4.5}
CHARGE = {'LYS':1.0,'ARG':1.0,'HIS':0.1,'ASP':-1.0,'GLU':-1.0}
HBD    = {'SER':1,'THR':1,'TYR':1,'TRP':1,'LYS':2,'ARG':3,'HIS':1,'ASN':1,'GLN':1,'CYS':1}
HBA    = {'ASP':2,'GLU':2,'SER':1,'THR':1,'TYR':1,'HIS':1,'ASN':1,'GLN':1,'MET':1}
AROM   = {'PHE','TYR','TRP','HIS'}
VDW_R  = {'C':1.70,'N':1.65,'O':1.60,'S':1.80,'P':1.80,
           'H':1.20,'F':1.47,'CL':1.75,'BR':1.85,'I':1.98}
WATER  = {'HOH','WAT','H2O'}

DESC_NAMES = ['hydrophobicity','charge','HBD','HBA','aromaticity',
              'linearity','planarity','sphericity','curvature']

def chem(atom):
    rn = atom.get_parent().get_resname().strip()
    return np.array([
        KD.get(rn, 0.0) / 4.5,
        CHARGE.get(rn, 0.0),
        HBD.get(rn, 0) / 3.0,
        HBA.get(rn, 0) / 2.0,
        1.0 if rn in AROM else 0.0,
    ], dtype=float)

def geom(coords, k=12):
    n = len(coords)
    shape, curv = np.zeros((n,3)), np.zeros(n)
    if n < 4:
        return shape, curv
    tree = KDTree(coords)
    kk = min(k, n-1)
    for i, p in enumerate(coords):
        _, idx = tree.query(p, k=kk+1)
        nb = coords[idx[1:]] - p
        if len(nb) < 3:
            continue
        cov = nb.T @ nb / len(nb)
        ev, _ = np.linalg.eigh(cov)
        lmin, lmid, lmax = ev[0], ev[1], ev[2]
        denom = lmax + 1e-12
        shape[i,0] = (lmax - lmid) / denom
        shape[i,1] = (lmid - lmin) / denom
        shape[i,2] = lmin          / denom
        curv[i]    = lmin / (ev.sum() + 1e-12)
    return shape, curv

def build_desc(atoms, coords):
    ch = np.array([chem(a) for a in atoms])
    sh, cv = geom(coords)
    return np.hstack([ch, sh, cv[:,None]])   # (N, 9)

def run_sasa(structure, n_pts=50):
    sr = ShrakeRupley(probe_radius=1.4, n_points=n_pts)
    sr.compute(structure, level='A')

def surface_atoms(structure, chain_ids=None, min_sasa=0.5):
    atoms, coords = [], []
    for model in structure:
        for chain in model:
            if chain_ids and chain.id not in chain_ids:
                continue
            for res in chain:
                if res.get_resname().strip() in WATER:
                    continue
                for atom in res:
                    if hasattr(atom, 'sasa') and atom.sasa > min_sasa:
                        atoms.append(atom)
                        coords.append(atom.coord.copy())
    arr = np.array(coords, dtype=float) if coords else np.empty((0,3))
    return atoms, arr

def get_chainA_ca(structure):
    """Cα coords + seq string for chain A (or first chain if no A)."""
    coords, seq = [], []
    for model in structure:
        chains = [c for c in model if c.id == 'A'] or list(model.get_chains())
        for chain in chains[:1]:
            for res in chain:
                if 'CA' in res and res.get_resname().strip() not in WATER:
                    coords.append(res['CA'].coord.copy())
                    seq.append(res.get_resname().strip()[:1])
        break
    return np.array(coords, dtype=np.float64), ''.join(seq)

# ── Reference setup ──────────────────────────────────────────────────────────

parser = PDBParser(QUIET=True)

print("Loading reference …")
ref_struct = parser.get_structure('ref', 'binding/A_B_b2_noPTM_A.pdb')

# Reference: SASA on chain A
ref_A = copy.deepcopy(ref_struct)
run_sasa(ref_A)
ref_atoms, ref_coords = surface_atoms(ref_A)
ref_desc = build_desc(ref_atoms, ref_coords)
ref_ca, ref_seq = get_chainA_ca(ref_struct)
print(f"  Reference surface atoms: {len(ref_atoms)}")

# Interface centroid from complex
cplx = parser.get_structure('cplx', 'binding/A_B_b2_noPTM.pdb')
A_ca_map, B_coords = {}, []
for model in cplx:
    for chain in model:
        for res in chain:
            if chain.id == 'A' and 'CA' in res:
                A_ca_map[(chain.id, res.get_id())] = res['CA'].coord.copy()
            elif chain.id == 'B':
                for atom in res: B_coords.append(atom.coord.copy())

tree_B = KDTree(np.array(B_coords))
iface_ca = np.array([c for k,c in A_ca_map.items()
                     if tree_B.query_ball_point(c, r=5.0)])
iface_centroid = iface_ca.mean(axis=0)
print(f"  Interface residues: {len(iface_ca)}  |  centroid: {np.round(iface_centroid,1)}")

IFACE_RADIUS = 10.0   # Å around centroid to select local surface atoms

# ── Groups ───────────────────────────────────────────────────────────────────

GROUP_TOP = {
    'A_B_b2_noPTM_A': 1,
    'A_B_b3_noPTM_A': 2,
    'b3_A':           3,
    '3d12_A':         4,
    'A_B_b3_A':       5,
    '6vy5_A':         6,
    '2vsm_A':         7,
    '6pdl_A':         8,
    '8ja5_A':         9,
    '2vsk_A':        10,
    '8xc4_A':        11,
    '3d11':          12,
    '9dvd_A':        13,
    'A_B_b2':        14,
    '6pd4_A':        15,
}

GROUP_G1 = {
    '2x9m_A': 24, '6cmg_A': 29, '8wjb_A': 32, '6p72_A': 33,
    '1v3d_A': 37, '7syy_A': 40, '8vf1_A': 42, '5nop_A': 44,
    '6vy4_A': 60, '2vwd_A': 63,
}

# ── Per-structure descriptor extraction ─────────────────────────────────────

def extract_iface_desc(name, rank, group_label):
    path = f'binding/{name}.pdb'
    struct = parser.get_structure(name, path)

    # For multi-chain structures, isolate chain A for SASA and descriptors
    st_A = copy.deepcopy(struct)
    for model in st_A:
        chains_to_remove = [c.id for c in model if c.id != 'A']
        for cid in chains_to_remove:
            model.detach_child(cid)
    # If no chain A exists, use the full structure as-is
    has_A = any(True for model in st_A for chain in model)
    target = st_A if has_A else struct

    run_sasa(target)
    atoms, coords = surface_atoms(target)
    if len(atoms) < 10:
        return None

    desc = build_desc(atoms, coords)

    # TM-align chain A Cα to reference
    cand_ca, cand_seq = get_chainA_ca(struct)
    if len(cand_ca) < 5:
        return None

    result   = tm_align(cand_ca, ref_ca, cand_seq, ref_seq)
    R, t     = result.u, result.t
    tm_score = result.tm_norm_chain2

    # Transform surface coords into reference frame
    coords_aln = (R @ coords.T).T + t

    # Select atoms near interface centroid in reference frame
    dists = np.linalg.norm(coords_aln - iface_centroid, axis=1)
    mask  = dists < IFACE_RADIUS
    if mask.sum() < 3:
        return None

    local_desc = desc[mask]   # (k, 9)
    return {
        'name':    name,
        'rank':    rank,
        'group':   group_label,
        'tm':      tm_score,
        'n_atoms': mask.sum(),
        'desc':    local_desc,                     # raw per-atom (k,9)
        'mean':    local_desc.mean(axis=0),        # (9,)
        'std':     local_desc.std(axis=0),         # (9,)
    }

print()
print("Processing TOP group (ranks 1–15) …")
top_records, top_all_desc = [], []
for name, rank in GROUP_TOP.items():
    r = extract_iface_desc(name, rank, 'TOP')
    if r:
        top_records.append(r)
        top_all_desc.append(r['desc'])
        print(f"  rank {rank:>2}  {name:<20}  TM={r['tm']:.4f}  n_iface_atoms={r['n_atoms']}")
    else:
        print(f"  rank {rank:>2}  {name:<20}  SKIPPED")

print()
print("Processing G1 group (ranks >16, same fold) …")
g1_records, g1_all_desc = [], []
for name, rank in GROUP_G1.items():
    r = extract_iface_desc(name, rank, 'G1')
    if r:
        g1_records.append(r)
        g1_all_desc.append(r['desc'])
        print(f"  rank {rank:>2}  {name:<20}  TM={r['tm']:.4f}  n_iface_atoms={r['n_atoms']}")
    else:
        print(f"  rank {rank:>2}  {name:<20}  SKIPPED")

# Pool all per-atom descriptors across each group
top_pool = np.vstack(top_all_desc)   # (N_top, 9)
g1_pool  = np.vstack(g1_all_desc)    # (N_g1, 9)

print(f"\n  TOP pooled atoms: {len(top_pool)}   G1 pooled atoms: {len(g1_pool)}")

# ── Per-structure mean summary ───────────────────────────────────────────────

print()
print("=" * 100)
print(f" Per-structure mean descriptors at interface region (within {IFACE_RADIUS} Å of interface centroid)")
print("=" * 100)
fmt_hdr = f"{'Grp':<4} {'Rank':>4}  {'Name':<20}  " + "  ".join(f"{n:>12}" for n in DESC_NAMES)
print(fmt_hdr)
print("-" * len(fmt_hdr))
for records in [top_records, g1_records]:
    for r in records:
        vals = "  ".join(f"{v:>12.4f}" for v in r['mean'])
        print(f"{r['group']:<4} {r['rank']:>4}  {r['name']:<20}  {vals}")
    print()

# ── Group-level comparison ───────────────────────────────────────────────────

print("=" * 100)
print(" Group comparison: TOP vs G1  (pooled per-atom descriptors)")
print("=" * 100)
print(f"{'Descriptor':<14}  {'TOP mean':>9} {'TOP std':>8}  {'G1 mean':>9} {'G1 std':>8}  "
      f"{'Δ(G1-TOP)':>10}  {'p-value':>9}  {'Verdict'}")
print("-" * 100)

rows_out = []
for i, dname in enumerate(DESC_NAMES):
    top_v = top_pool[:, i]
    g1_v  = g1_pool[:, i]
    delta = g1_v.mean() - top_v.mean()
    stat, pval = mannwhitneyu(top_v, g1_v, alternative='two-sided')
    sig = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else 'ns'))
    print(f"{dname:<14}  {top_v.mean():>9.4f} {top_v.std():>8.4f}  "
          f"{g1_v.mean():>9.4f} {g1_v.std():>8.4f}  "
          f"{delta:>+10.4f}  {pval:>9.2e}  {sig}")
    rows_out.append(dict(descriptor=dname,
                         top_mean=top_v.mean(), top_std=top_v.std(),
                         g1_mean=g1_v.mean(),  g1_std=g1_v.std(),
                         delta=delta, pval=pval, sig=sig))

# ── Save TSV ─────────────────────────────────────────────────────────────────

tsv_out = 'results/dmasif/fp_comparison.tsv'
with open(tsv_out, 'w') as f:
    f.write('descriptor\ttop_mean\ttop_std\tg1_mean\tg1_std\tdelta_g1_minus_top\tpval\tsig\n')
    for r in rows_out:
        f.write(f"{r['descriptor']}\t{r['top_mean']:.5f}\t{r['top_std']:.5f}\t"
                f"{r['g1_mean']:.5f}\t{r['g1_std']:.5f}\t"
                f"{r['delta']:+.5f}\t{r['pval']:.3e}\t{r['sig']}\n")

# Also save per-structure means
tsv_per = 'results/dmasif/fp_per_structure.tsv'
with open(tsv_per, 'w') as f:
    f.write('group\trank\tname\ttm_score\tn_iface_atoms\t' + '\t'.join(DESC_NAMES) + '\n')
    for records in [top_records, g1_records]:
        for r in records:
            f.write(f"{r['group']}\t{r['rank']}\t{r['name']}\t{r['tm']:.4f}\t{r['n_atoms']}\t"
                    + '\t'.join(f"{v:.5f}" for v in r['mean']) + '\n')

print(f"\nSaved → {tsv_out}")
print(f"Saved → {tsv_per}")
print()
print("Legend:")
print(f"  Δ = G1 mean − TOP mean  (positive = G1 higher, negative = G1 lower)")
print(f"  *** p<0.001  ** p<0.01  * p<0.05  ns = not significant")
print(f"  Interface window: {IFACE_RADIUS} Å radius around reference interface centroid")
