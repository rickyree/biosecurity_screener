#!/usr/bin/env python3
"""
dMaSIF-inspired protein surface interaction fingerprinting.

Workflow:
  1. Load A-B complex (binding/A_B_b2.pdb), compute whole-complex SASA.
  2. Identify chain-A surface atoms at the A-B interface (<= 5 Å from chain B).
  3. Compute geometric + chemical descriptors and 18-dim patch fingerprints
     for each interface atom on A.
  4. For each candidate C protein:
       a. Compute outer-surface atoms (SASA on the whole multi-chain structure).
       b. Build patch fingerprints for every C surface atom.
       c. Match each A-interface fingerprint to its nearest C fingerprint (NN).
       d. RANSAC rigid-body alignment: src = C matched coords, tgt = A interface coords.
       e. Apply transform to all C atoms → docked pose.
       f. Steric clash check against chain B.
  5. Report ranked predictions.
"""

import os, sys, copy, time, warnings, argparse
import numpy as np
from scipy.spatial import KDTree
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings('ignore')
np.random.seed(42)

os.makedirs('results/dmasif', exist_ok=True)

A_CACHE        = 'results/dmasif/A_interface_cache.npz'
A_CACHE_ESM    = 'results/dmasif/A_interface_cache_esmfold.npz'
A_CACHE_APO    = 'results/dmasif/A_interface_cache_apo_nosas.npz'
A_CACHE_SGK    = 'results/dmasif/A_interface_cache_sgk.npz'
ESMFOLD_A_PDB  = 'results/esmfold/A_standalone/esmfold_esmfold_A_standalone_structure.pdb'
APO_A_PDB      = 'binding/2vwd_A.pdb'
SGK_A_PDB      = 'binding/1sgk_A.pdb'   # apo receptor for NAD-site screen
F0L_A_PDB      = 'binding/1f0l_A.pdb'   # holo reference (chain A + APU)
F0L_CHAIN      = 'A'
F0L_LIGAND     = 'APU'
PLDDT_THRESH   = 50.0   # minimum ESMFold pLDDT to include an interface residue
EXCLUDE_RNUMS  = {579, 580}  # induced-fit loop residues (ordered only in bound state)

# ── Residue property tables ──────────────────────────────────────────────────

# Kyte-Doolittle hydrophobicity, normalised to [-1, 1] (divide by 4.5)
KD = {
    'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
    'MET': 1.9, 'ALA': 1.8, 'GLY':-0.4, 'THR':-0.7, 'SER':-0.8,
    'TRP':-0.9, 'TYR':-1.3, 'PRO':-1.6, 'HIS':-3.2, 'GLU':-3.5,
    'GLN':-3.5, 'ASP':-3.5, 'ASN':-3.5, 'LYS':-3.9, 'ARG':-4.5,
}
CHARGE  = {'LYS': 1.0, 'ARG': 1.0, 'HIS': 0.1, 'ASP':-1.0, 'GLU':-1.0}
HBD     = {'SER':1,'THR':1,'TYR':1,'TRP':1,'LYS':2,'ARG':3,'HIS':1,'ASN':1,'GLN':1,'CYS':1}
HBA     = {'ASP':2,'GLU':2,'SER':1,'THR':1,'TYR':1,'HIS':1,'ASN':1,'GLN':1,'MET':1}
AROM    = {'PHE','TYR','TRP','HIS'}
VDW_R   = {'C':1.70,'N':1.65,'O':1.60,'S':1.80,'P':1.80,
            'H':1.20,'F':1.47,'CL':1.75,'BR':1.85,'I':1.98}
WATER   = {'HOH','WAT','H2O'}

def vdw(atom):
    elem = (atom.element or atom.name[0]).strip().upper()
    return VDW_R.get(elem, 1.70)

# ── Surface extraction ───────────────────────────────────────────────────────

def run_sasa(structure, n_pts=50):
    sr = ShrakeRupley(probe_radius=1.4, n_points=n_pts)
    sr.compute(structure, level='A')

def surface_atoms(structure, chain_ids=None, min_sasa=0.5):
    """Return (atoms[], coords ndarray) for SASA-exposed non-water atoms."""
    atoms, coords = [], []
    for model in structure:
        for chain in model:
            if chain_ids and chain.id not in chain_ids:
                continue
            for res in chain:
                rn = res.get_resname().strip()
                if rn in WATER:
                    continue
                for atom in res:
                    if hasattr(atom, 'sasa') and atom.sasa > min_sasa:
                        atoms.append(atom)
                        coords.append(atom.coord.copy())
    arr = np.array(coords, dtype=float) if coords else np.empty((0, 3))
    return atoms, arr

def all_heavy_atoms(structure, chain_ids=None):
    """Return (atoms[], coords ndarray) for all non-water, non-hydrogen atoms."""
    atoms, coords = [], []
    for model in structure:
        for chain in model:
            if chain_ids and chain.id not in chain_ids:
                continue
            for res in chain:
                if res.get_resname().strip() in WATER:
                    continue
                for atom in res:
                    elem = (atom.element or atom.name[0]).strip().upper()
                    if elem == 'H':
                        continue
                    atoms.append(atom)
                    coords.append(atom.coord.copy())
    arr = np.array(coords, dtype=float) if coords else np.empty((0, 3))
    return atoms, arr

# ── Per-atom chemical features (5-dim) ──────────────────────────────────────

def chem(atom):
    rn = atom.get_parent().get_resname().strip()
    return np.array([
        KD.get(rn, 0.0) / 4.5,          # hydrophobicity   [-1,1]
        CHARGE.get(rn, 0.0),             # formal charge    {-1,0,+1}
        HBD.get(rn, 0)    / 3.0,         # H-bond donors    [0,1]
        HBA.get(rn, 0)    / 2.0,         # H-bond acceptors [0,1]
        1.0 if rn in AROM else 0.0,      # aromaticity      {0,1}
    ], dtype=float)

# ── Geometric features: normals + curvature via local PCA ───────────────────

def geom(coords, k=12):
    """Return rotation-invariant shape descriptors (N,3) and curvature (N,).

    Replaces raw normal vectors (which are frame-dependent) with eigenvalue
    ratios from local PCA — these describe surface shape regardless of the
    global orientation of the structure:
      shape[:,0] = linearity  = (λmax - λmid) / (λmax + ε)  — ridge-like
      shape[:,1] = planarity  = (λmid - λmin) / (λmax + ε)  — flat surface
      shape[:,2] = sphericity = λmin          / (λmax + ε)  — dome/bump
    curvature   = λmin / (Σλ + ε)  — fraction of variance along normal
    All four are invariant to rotation and translation.
    """
    n = len(coords)
    shape     = np.zeros((n, 3))
    curvature = np.zeros(n)
    if n < 4:
        return shape, curvature
    tree = KDTree(coords)
    kk   = min(k, n - 1)
    for i, p in enumerate(coords):
        _, idx = tree.query(p, k=kk + 1)
        nb = coords[idx[1:]] - p
        if len(nb) < 3:
            continue
        cov = nb.T @ nb / len(nb)
        ev, _ = np.linalg.eigh(cov)          # ascending: ev[0] ≤ ev[1] ≤ ev[2]
        lmin, lmid, lmax = ev[0], ev[1], ev[2]
        denom = lmax + 1e-12
        shape[i, 0] = (lmax - lmid) / denom  # linearity
        shape[i, 1] = (lmid - lmin) / denom  # planarity
        shape[i, 2] = lmin          / denom  # sphericity
        curvature[i] = lmin / (ev.sum() + 1e-12)
    return shape, curvature

# ── Full descriptor matrix (N, 9) ────────────────────────────────────────────

def build_desc(atoms, coords):
    ch       = np.array([chem(a) for a in atoms])  # (N, 5)
    sh, cv   = geom(coords)                         # (N, 3),  (N,)
    return np.hstack([ch, sh, cv[:, None]])          # (N, 9)

# ── Patch fingerprint (18-dim) ───────────────────────────────────────────────

def patch_fp(idx, coords, desc, tree, radius=12.0):
    """Weighted-mean + std of descriptors within `radius` Å."""
    p    = coords[idx]
    near = tree.query_ball_point(p, r=radius)
    if len(near) < 3:
        return None
    near  = np.array(near)
    d     = np.linalg.norm(coords[near] - p, axis=1)
    w     = 1.0 - d / radius
    w    /= w.sum()
    mu    = (desc[near] * w[:, None]).sum(0)
    sd    = desc[near].std(0)
    return np.concatenate([mu, sd])              # 18-dim

def build_all_fps(atoms, coords, desc, radius=12.0, max_pts=3000):
    """Compute patch fingerprints for all (or subsampled) surface atoms.
    Returns fps (M, 18), fp_to_atom_idx (M,).
    """
    n    = len(atoms)
    idxs = np.arange(n)

    # Subsample very large surfaces for speed
    if n > max_pts:
        step  = n // max_pts
        idxs  = idxs[::step]

    tree = KDTree(coords)
    fps, valid = [], []
    for i in idxs:
        fp = patch_fp(i, coords, desc, tree, radius)
        if fp is not None:
            fps.append(fp)
            valid.append(i)
    return np.array(fps), np.array(valid)

# ── Kabsch rigid-body alignment ──────────────────────────────────────────────

def kabsch(P, Q):
    """R, t  s.t.  R @ P[i] + t  ≈  Q[i]  (least squares)."""
    cp, cq   = P.mean(0), Q.mean(0)
    Pc, Qc   = P - cp,   Q - cq
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d        = np.linalg.det(Vt.T @ U.T)
    R        = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp

def ransac_rigid(src, tgt, n_iter=600, thr=3.0, min_k=3):
    """RANSAC rigid-body: src[i] → tgt[i] correspondences.
    Returns R, t, inlier_bool_array.
    """
    m = len(src)
    if m < min_k:
        return np.eye(3), np.zeros(3), np.zeros(m, bool)
    best_R, best_t = np.eye(3), np.zeros(3)
    best_mask = np.zeros(m, bool)
    best_n    = 0
    for _ in range(n_iter):
        s = np.random.choice(m, min_k, replace=False)
        try:
            R, t = kabsch(src[s], tgt[s])
        except Exception:
            continue
        dists = np.linalg.norm((R @ src.T).T + t - tgt, axis=1)
        mask  = dists < thr
        if mask.sum() > best_n:
            best_n    = mask.sum()
            best_mask = mask.copy()
            if mask.sum() >= min_k:
                try:    best_R, best_t = kabsch(src[mask], tgt[mask])
                except: best_R, best_t = R, t
            else:
                best_R, best_t = R, t
    return best_R, best_t, best_mask

# ── Steric clash check ───────────────────────────────────────────────────────

def clash_check(c_atoms, c_pos, b_pos, b_vdw_arr, factor=0.6):
    """Count C atoms that clash with any B atom (VdW overlap > factor).
    b_vdw_arr: precomputed VdW radii array for B surface atoms (N_B,).
    """
    tree  = KDTree(b_pos)
    clash = 0
    for a, p in zip(c_atoms, c_pos):
        ra   = vdw(a)
        near = tree.query_ball_point(p, r=ra + 2.0)
        for j in near:
            if np.linalg.norm(p - b_pos[j]) < factor * (ra + b_vdw_arr[j]):
                clash += 1
                break
    return clash, clash / max(len(c_atoms), 1)

# ── Save transformed structure ───────────────────────────────────────────────

def save_docked(structure, R, t, path):
    s = copy.deepcopy(structure)
    for a in s.get_atoms():
        a.set_coord(R @ a.coord + t)
    io = PDBIO()
    io.set_structure(s)
    io.save(path)

# ── Main ─────────────────────────────────────────────────────────────────────

def superimpose_kabsch(mobile_ca, target_ca):
    """Return R, t that maps mobile CA coords onto target CA coords (Kabsch)."""
    R, t = kabsch(mobile_ca, target_ca)
    return R, t


def load_or_build_A_fingerprints(pdb_parser, use_esmfold=False, use_apo=False, use_sgk=False):
    """Return (A_fps, A_iface_coords, B_pos, B_vdw) either from cache or by computing.

    use_sgk:     fingerprints from 1sgk_A (apo DT) at NAD-site residues identified
                 by superimposing onto 1f0l_A and projecting APU coordinates.
                 B_pos/B_vdw are empty (no partner — clash check is skipped).
    use_esmfold: fingerprints from ESMFold standalone A (pLDDT-filtered).
    use_apo:     fingerprints from apo crystal monomer (2vwd_A), superimposed
                 onto bound frame; excludes induced-fit residues (EXCLUDE_RNUMS).
    default:     fingerprints from bound complex chain A.
    """
    cache = A_CACHE_SGK if use_sgk else (A_CACHE_APO if use_apo else (A_CACHE_ESM if use_esmfold else A_CACHE))
    if os.path.exists(cache):
        print(f"  Loading cached A-interface fingerprints from {cache}")
        d = np.load(cache)
        return d['A_fps'], d['A_iface_coords'], d['B_pos'], d['B_vdw']

    print("=" * 65)
    print(" Step 1: Building dMaSIF fingerprint for A's interface with B")
    if use_sgk:
        print(f"         (1sgk_A apo, NAD-site via superposition onto 1f0l_A)")
    elif use_esmfold:
        print(f"         (ESMFold standalone, pLDDT>{PLDDT_THRESH:.0f} filter)")
    elif use_apo:
        print(f"         (apo crystal monomer 2vwd_A, exclude residues {sorted(EXCLUDE_RNUMS)})")
    print("=" * 65)

    if use_sgk:
        # ── 1sgk_A NAD-site path ──────────────────────────────────────────────
        sgk_struct   = pdb_parser.get_structure('SGK', SGK_A_PDB)
        f0l_struct   = pdb_parser.get_structure('F0L', F0L_A_PDB)

        # Cα coords and sequences of 1f0l chain A (protein residues)
        from Bio.SeqUtils import seq1
        from Bio import pairwise2
        f0l_residues = [res for model in f0l_struct for chain in model
                        if chain.id == F0L_CHAIN
                        for res in chain if res.id[0] == ' ' and 'CA' in res]
        f0l_ca  = np.array([res['CA'].coord for res in f0l_residues])
        f0l_seq = ''.join(seq1(r.get_resname()) for r in f0l_residues)

        # Cα coords and sequence of 1sgk_A
        sgk_chain_id = list(sgk_struct[0].get_chains())[0].id
        sgk_residues = [res for model in sgk_struct for chain in model
                        for res in chain if res.id[0] == ' ' and 'CA' in res]
        sgk_ca  = np.array([res['CA'].coord for res in sgk_residues])
        sgk_seq = ''.join(seq1(r.get_resname()) for r in sgk_residues)

        # Sequence alignment → matched Cα pairs for Kabsch
        aln = pairwise2.align.globalds(
            sgk_seq, f0l_seq,
            pairwise2.substitution_matrices.load("BLOSUM62"),
            -10, -0.5, one_alignment_only=True)[0]
        sgk_aln, f0l_aln = aln.seqA, aln.seqB
        sgk_idx, f0l_idx = [], []
        si, fi = 0, 0
        for a, b in zip(sgk_aln, f0l_aln):
            if a != '-' and b != '-':
                sgk_idx.append(si)
                f0l_idx.append(fi)
            if a != '-': si += 1
            if b != '-': fi += 1
        sgk_ca_aln = sgk_ca[sgk_idx]
        f0l_ca_aln = f0l_ca[f0l_idx]

        # Kabsch: superimpose 1sgk_A onto 1f0l chain A using sequence-aligned pairs
        R_sup, t_sup = superimpose_kabsch(sgk_ca_aln, f0l_ca_aln)
        rmsd_sup = float(np.sqrt(np.mean(
            np.sum(((R_sup @ sgk_ca_aln.T).T + t_sup - f0l_ca_aln)**2, axis=1))))
        print(f"  Superimposed 1sgk_A onto 1f0l chain A "
              f"({len(sgk_idx)} aligned Cα pairs, RMSD={rmsd_sup:.2f} Å)")

        # APU atoms in 1f0l chain A → inverse-transform into 1sgk_A frame
        apu_coords = np.array([a.coord
                                for model in f0l_struct for chain in model
                                if chain.id == F0L_CHAIN
                                for res in chain for a in res
                                if res.get_resname().strip() == F0L_LIGAND])
        print(f"  APU atoms in 1f0l: {len(apu_coords)}")
        R_inv = R_sup.T
        t_inv = -R_inv @ t_sup
        apu_in_sgk = (R_inv @ apu_coords.T).T + t_inv

        apu_tree = KDTree(apu_in_sgk)
        iface_res = []
        for res in sgk_residues:
            for atom in res:
                if apu_tree.query_ball_point(atom.coord, r=5.0):
                    iface_res.append(res)
                    break
        print(f"  1sgk_A interface residues ({len(iface_res)}): {[r.id[1] for r in iface_res]}")

        A_atoms, A_pos = all_heavy_atoms(sgk_struct)
        print(f"  1sgk_A heavy atoms: {len(A_atoms)}")
        A_desc = build_desc(A_atoms, A_pos)
        A_tree = KDTree(A_pos)

        A_fps, A_fp_coords = [], []
        for res in iface_res:
            if 'CA' not in res:
                continue
            ca_coord = res['CA'].coord
            _, idx = A_tree.query(ca_coord)
            fp = patch_fp(idx, A_pos, A_desc, A_tree)
            if fp is not None:
                A_fps.append(fp)
                A_fp_coords.append(A_pos[idx])

        # No partner B — return empty arrays; clash check will be skipped
        B_pos = np.empty((0, 3))
        B_vdw = np.empty(0)

        if len(A_fps) == 0:
            sys.exit("ERROR: no valid fingerprints computed for 1sgk_A.")

        A_fps          = np.array(A_fps)
        A_iface_coords = np.array(A_fp_coords)
        print(f"  Interface fingerprints: {len(A_fps)}")

        np.savez(cache, A_fps=A_fps, A_iface_coords=A_iface_coords,
                 B_pos=B_pos, B_vdw=B_vdw)
        print(f"  Cached A fingerprints → {cache}")
        return A_fps, A_iface_coords, B_pos, B_vdw

    ab = pdb_parser.get_structure('AB', 'binding/A_B_b2_noPTM.pdb')

    # ── B surface (needed for interface detection and clash check) ────────────
    print("  Computing SASA for chain B (partner) …", end=' ', flush=True)
    ab_B_only = copy.deepcopy(ab)
    for model in ab_B_only:
        for cid in [c.id for c in model if c.id != 'B']:
            model.detach_child(cid)
    run_sasa(ab_B_only)
    print("done")
    B_sasa_map = {}
    for model in ab_B_only:
        for chain in model:
            for res in chain:
                for atom in res:
                    B_sasa_map[(res.get_full_id(), atom.name)] = atom.sasa
    for model in ab:
        for chain in model:
            if chain.id != 'B':
                continue
            for res in chain:
                for atom in res:
                    key = (res.get_full_id(), atom.name)
                    if key in B_sasa_map:
                        atom.sasa = B_sasa_map[key]

    B_atoms, B_pos = surface_atoms(ab, chain_ids=['B'])
    print(f"  Chain B surface atoms : {len(B_atoms)}")
    B_vdw = np.array([vdw(a) for a in B_atoms])

    # ── Identify interface residue positions (sequential, 0-based) ───────────
    # Done on the bound complex so we know which positions contact B.
    tree_B_all = KDTree(np.array([a.coord for model in ab
                                  for chain in model if chain.id == 'B'
                                  for res in chain for a in res]))
    ab_a_residues = [res for model in ab for chain in model if chain.id == 'A'
                     for res in chain if res.id[0] == ' ']
    iface_seq_pos = set()
    for i, res in enumerate(ab_a_residues):
        for atom in res:
            if tree_B_all.query_ball_point(atom.coord, r=5.0):
                iface_seq_pos.add(i)
                break

    # Map sequential position → original residue number (for exclusion)
    pos_to_rnum = {i: res.id[1] for i, res in enumerate(ab_a_residues)}

    if use_apo:
        # ── Apo crystal monomer path ──────────────────────────────────────────
        if not os.path.exists(APO_A_PDB):
            sys.exit(f"ERROR: Apo structure not found at {APO_A_PDB}")

        apo_struct   = pdb_parser.get_structure('APO', APO_A_PDB)
        apo_residues = [res for model in apo_struct for chain in model
                        for res in chain if res.id[0] == ' ' and 'CA' in res]

        # Exclude induced-fit residues by original residue number
        good_pos = [p for p in sorted(iface_seq_pos)
                    if pos_to_rnum.get(p) not in EXCLUDE_RNUMS]
        print(f"  Interface positions (bound complex) : {len(iface_seq_pos)}")
        print(f"  After excluding residues {sorted(EXCLUDE_RNUMS)}: {len(good_pos)}")

        # Superimpose apo onto bound chain A
        n_shared = min(len(apo_residues), len(ab_a_residues))
        apo_ca = np.array([apo_residues[i]['CA'].coord for i in range(n_shared)])
        ab_ca  = np.array([ab_a_residues[i]['CA'].coord for i in range(n_shared)
                           if 'CA' in ab_a_residues[i]])
        n_aln  = min(len(apo_ca), len(ab_ca))
        R_sup, t_sup = superimpose_kabsch(apo_ca[:n_aln], ab_ca[:n_aln])
        print(f"  Superimposed apo-A onto bound chain A ({n_aln} Cα pairs)")

        apo_sup = copy.deepcopy(apo_struct)
        for atom in apo_sup.get_atoms():
            atom.set_coord(R_sup @ atom.coord + t_sup)

        A_atoms, A_pos = all_heavy_atoms(apo_sup)
        print(f"  Apo-A heavy atoms : {len(A_atoms)}")

        A_desc = build_desc(A_atoms, A_pos)
        A_tree = KDTree(A_pos)
        apo_sup_residues = [res for model in apo_sup for chain in model
                            for res in chain if res.id[0] == ' ' and 'CA' in res]
        A_fps, A_fp_coords = [], []
        for p in good_pos:
            if p >= len(apo_sup_residues):
                continue
            ca_coord = apo_sup_residues[p]['CA'].coord
            # nearest heavy atom to Cα is the Cα itself — no SASA filter needed
            _, idx = A_tree.query(ca_coord)
            fp = patch_fp(idx, A_pos, A_desc, A_tree)
            if fp is not None:
                A_fps.append(fp)
                A_fp_coords.append(A_pos[idx])

    elif use_esmfold:
        # ── ESMFold path ──────────────────────────────────────────────────────
        if not os.path.exists(ESMFOLD_A_PDB):
            sys.exit(f"ERROR: ESMFold structure not found at {ESMFOLD_A_PDB}")

        esm_struct = pdb_parser.get_structure('ESM', ESMFOLD_A_PDB)
        esm_residues = [res for model in esm_struct for chain in model
                        for res in chain if res.id[0] == ' ' and 'CA' in res]

        # Filter interface positions to well-predicted residues only
        good_pos = [p for p in sorted(iface_seq_pos)
                    if p < len(esm_residues)
                    and esm_residues[p]['CA'].get_bfactor() >= PLDDT_THRESH]
        print(f"  Interface positions (bound complex) : {len(iface_seq_pos)}")
        print(f"  After ESMFold pLDDT>{PLDDT_THRESH:.0f} filter        : {len(good_pos)}")

        # Superimpose ESMFold-A onto bound chain A via shared CA positions
        n_shared = min(len(esm_residues), len(ab_a_residues))
        esm_ca  = np.array([esm_residues[i]['CA'].coord for i in range(n_shared)])
        ab_ca   = np.array([ab_a_residues[i]['CA'].coord
                            for i in range(n_shared) if 'CA' in ab_a_residues[i]])
        # Both must have same length for Kabsch
        n_aln = min(len(esm_ca), len(ab_ca))
        R_sup, t_sup = superimpose_kabsch(esm_ca[:n_aln], ab_ca[:n_aln])
        print(f"  Superimposed ESMFold-A onto bound chain A ({n_aln} Cα pairs)")

        # Apply superposition to all ESMFold atoms
        esm_sup = copy.deepcopy(esm_struct)
        for atom in esm_sup.get_atoms():
            atom.set_coord(R_sup @ atom.coord + t_sup)

        # Compute SASA on superimposed ESMFold-A (standalone)
        print("  Computing SASA for ESMFold-A …", end=' ', flush=True)
        run_sasa(esm_sup)
        print("done")

        A_atoms, A_pos = surface_atoms(esm_sup)
        print(f"  ESMFold-A surface atoms : {len(A_atoms)}")

        # Build descriptors on the full ESMFold surface
        A_desc = build_desc(A_atoms, A_pos)
        A_tree = KDTree(A_pos)

        # Fingerprint only the good interface positions: find nearest ESMFold
        # surface atom to each good-position CA (after superposition)
        esm_sup_residues = [res for model in esm_sup for chain in model
                            for res in chain if res.id[0] == ' ' and 'CA' in res]
        surf_tree = KDTree(A_pos)
        A_fps, A_fp_coords = [], []
        for p in good_pos:
            if p >= len(esm_sup_residues):
                continue
            ca_coord = esm_sup_residues[p]['CA'].coord
            dist, idx = surf_tree.query(ca_coord)
            if dist > 5.0:          # CA not on surface — skip
                continue
            fp = patch_fp(idx, A_pos, A_desc, A_tree)
            if fp is not None:
                A_fps.append(fp)
                A_fp_coords.append(A_pos[idx])

    else:
        # ── Bound complex path (original behaviour) ───────────────────────────
        ab_A_only = copy.deepcopy(ab)
        for model in ab_A_only:
            for cid in [c.id for c in model if c.id != 'A']:
                model.detach_child(cid)
        print("  Computing SASA for chain A alone …", end=' ', flush=True)
        run_sasa(ab_A_only)
        print("done")

        sasa_map = {}
        for model in ab_A_only:
            for chain in model:
                for res in chain:
                    for atom in res:
                        sasa_map[(res.get_full_id(), atom.name)] = atom.sasa
        for model in ab:
            for chain in model:
                if chain.id != 'A':
                    continue
                for res in chain:
                    for atom in res:
                        key = (res.get_full_id(), atom.name)
                        if key in sasa_map:
                            atom.sasa = sasa_map[key]

        A_atoms, A_pos = surface_atoms(ab, chain_ids=['A'])
        print(f"  Chain A surface atoms : {len(A_atoms)}  (chain-A-alone SASA)")

        if len(A_atoms) == 0:
            sys.exit("ERROR: could not extract surface atoms.")

        tree_B     = KDTree(B_pos)
        iface_mask = np.array([bool(tree_B.query_ball_point(p, r=5.0)) for p in A_pos])
        n_iface    = iface_mask.sum()
        print(f"  Interface atoms on A  : {n_iface}  (≤5 Å from chain B)")
        if n_iface == 0:
            sys.exit("ERROR: no interface atoms at 5 Å cutoff.")

        A_desc = build_desc(A_atoms, A_pos)
        A_tree = KDTree(A_pos)
        A_fps, A_fp_coords = [], []
        for i in np.where(iface_mask)[0]:
            fp = patch_fp(i, A_pos, A_desc, A_tree)
            if fp is not None:
                A_fps.append(fp)
                A_fp_coords.append(A_pos[i])

    if len(A_fps) == 0:
        sys.exit("ERROR: no valid fingerprints computed.")

    A_fps          = np.array(A_fps)
    A_iface_coords = np.array(A_fp_coords)
    print(f"  Interface fingerprints : {len(A_fps)}")

    np.savez(cache,
             A_fps=A_fps, A_iface_coords=A_iface_coords,
             B_pos=B_pos, B_vdw=B_vdw)
    print(f"  Cached A fingerprints → {cache}")

    return A_fps, A_iface_coords, B_pos, B_vdw


def screen_candidates(candidates, A_fps, A_iface_coords, B_pos, B_vdw, pdb_parser):
    results = []
    for c_path in candidates:
        name = os.path.basename(c_path).replace('.pdb', '')
        print()
        print("─" * 65)
        print(f" Candidate C : {name}")
        print("─" * 65)
        tc = time.time()

        c_struct = pdb_parser.get_structure('C', c_path)

        C_atoms, C_pos = all_heavy_atoms(c_struct)
        print(f"  C heavy atoms : {len(C_atoms)}")
        if len(C_atoms) < 10:
            print("  SKIP: too few atoms")
            continue

        print("  Building descriptors for C …", end=' ', flush=True)
        C_desc = build_desc(C_atoms, C_pos)
        C_tree = KDTree(C_pos)
        print("done")

        print("  Building Cα-anchored patch fingerprints for C …", end=' ', flush=True)
        C_fps_list, C_fp_idx_list = [], []
        for i, atom in enumerate(C_atoms):
            if atom.name == 'CA':
                fp = patch_fp(i, C_pos, C_desc, C_tree)
                if fp is not None:
                    C_fps_list.append(fp)
                    C_fp_idx_list.append(i)
        C_fps    = np.array(C_fps_list)    if C_fps_list    else np.empty((0, 18))
        C_fp_idx = np.array(C_fp_idx_list) if C_fp_idx_list else np.empty(0, int)
        print(f"{len(C_fps)} fingerprints")
        if len(C_fps) == 0:
            print("  SKIP: no valid fingerprints")
            continue

        # NN matching in fingerprint space
        print("  NN matching A-interface FPs → C surface …", end=' ', flush=True)
        tree_cfp         = KDTree(C_fps)
        nn_dists, nn_idx = tree_cfp.query(A_fps)
        top_k            = max(3, len(nn_dists) // 3)
        top_sel          = np.argsort(nn_dists)[:top_k]
        mean_fp_dist     = nn_dists[top_sel].mean()
        print(f"done  (mean top-30% FP dist = {mean_fp_dist:.4f})")

        src_corr = C_pos[C_fp_idx[nn_idx[top_sel]]]
        tgt_corr = A_iface_coords[top_sel]

        # RANSAC rigid alignment
        print("  RANSAC rigid alignment …", end=' ', flush=True)
        R, t, inlier_mask = ransac_rigid(src_corr, tgt_corr, n_iter=600, thr=3.0)
        n_in   = inlier_mask.sum()
        in_rat = n_in / max(len(inlier_mask), 1)
        if n_in >= 3:
            rmsd = float(np.sqrt(np.mean(
                np.sum(((R @ src_corr[inlier_mask].T).T + t - tgt_corr[inlier_mask])**2, axis=1)
            )))
        else:
            rmsd = float('inf')
        print(f"done  (inliers {n_in}/{len(inlier_mask)}, patch RMSD {rmsd:.2f} Å)")

        # Dock + clash check (skipped when no partner B)
        all_C_atoms = list(c_struct.get_atoms())
        all_C_pos   = np.array([a.coord for a in all_C_atoms])
        C_docked    = (R @ all_C_pos.T).T + t

        if len(B_pos) > 0:
            print("  Steric clash check …", end=' ', flush=True)
            n_clash, clash_r = clash_check(all_C_atoms, C_docked, B_pos, B_vdw)
            print(f"done  ({n_clash} clashes, {clash_r:.1%} of C atoms)")
        else:
            n_clash, clash_r = 0, 0.0

        out_pdb = f"results/dmasif/{name}_docked.pdb"
        save_docked(c_struct, R, t, out_pdb)
        print(f"  Saved docked pose → {out_pdb}  ({time.time()-tc:.0f}s)")

        results.append({
            'name'      : name,
            'fp_dist'   : mean_fp_dist,
            'inlier_n'  : n_in,
            'inlier_rat': in_rat,
            'rmsd'      : rmsd,
            'n_clash'   : n_clash,
            'clash_r'   : clash_r,
        })
    return results


def print_summary(results, tsv_path='results/dmasif/rankings.tsv'):
    if not results:
        return

    print()
    print("=" * 72)
    print(" RESULTS — dMaSIF-based C-B interaction screen")
    print("=" * 72)
    print(f"{'Candidate':<26} {'FP-dist':>8} {'Inlier%':>8} {'RMSD(Å)':>8} {'Clashes':>8} {'Clash%':>7}")
    print("-" * 72)
    for r in sorted(results, key=lambda x: x['fp_dist']):
        rmsd_str = f"{r['rmsd']:8.2f}" if np.isfinite(r['rmsd']) else "      — "
        print(f"{r['name']:<26} {r['fp_dist']:>8.4f} {r['inlier_rat']:>8.1%} "
              f"{rmsd_str} {r['n_clash']:>8d} {r['clash_r']:>7.1%}")
    print()
    print("Legend:")
    print("  FP-dist  lower → better surface match to A's interface fingerprint")
    print("  Inlier%  higher → more consistent RANSAC 3-D alignment")
    print("  RMSD     lower  → tighter patch fit onto A's interface")
    print("  Clash%   lower  → cleaner dock against B (no steric overlaps)")

    def norm01(arr):
        lo, hi = arr.min(), arr.max()
        return (arr - lo) / (hi - lo + 1e-12)

    fp_arr = np.array([r['fp_dist']    for r in results])
    cr_arr = np.array([r['clash_r']    for r in results])
    ir_arr = np.array([r['inlier_rat'] for r in results])
    raw_rd = np.array([r['rmsd']       for r in results])
    finite = raw_rd[np.isfinite(raw_rd)]
    penalty = finite.max() * 2.0 if len(finite) > 0 else 10.0
    rd_arr  = np.where(np.isfinite(raw_rd), raw_rd, penalty)

    # Score: surface similarity + alignment quality (no steric clash term)
    score = norm01(fp_arr) + norm01(rd_arr) - norm01(ir_arr)

    ranked = sorted(zip(score, results), key=lambda x: x[0])
    median = np.median(score)

    W = 30
    print()
    print("=" * 72)
    print(" RANKING — surface fingerprint match + RANSAC alignment quality")
    print("=" * 72)
    print(f"  {'Rank':<6} {'Candidate':<{W}} {'Score':>7}  {'FP-dist':>8}  {'Inlier%':>8}  {'RMSD':>6}")
    print(f"  {'-'*6} {'-'*W} {'-'*7}  {'-'*8}  {'-'*8}  {'-'*6}")
    pass  # tsv_path comes from parameter
    with open(tsv_path, 'w') as tsv:
        tsv.write('rank\tcandidate\tscore\tfp_dist\tinlier_pct\trmsd_A\tn_clash\tclash_pct\tverdict\n')
        for rank, (sc, r) in enumerate(ranked, 1):
            rmsd_s = f"{r['rmsd']:.2f}" if np.isfinite(r['rmsd']) else "—"
            verdict = "PREDICTED BINDER" if sc <= median else "unlikely binder"
            print(f"  {rank:<6} {r['name']:<{W}} {sc:>+7.3f}  {r['fp_dist']:>8.4f}  {r['inlier_rat']:>8.1%}  {rmsd_s:>6}   [{verdict}]")
            tsv.write(f"{rank}\t{r['name']}\t{sc:+.4f}\t{r['fp_dist']:.4f}\t{r['inlier_rat']:.4f}\t{rmsd_s}\t{r['n_clash']}\t{r['clash_r']:.4f}\t{verdict}\n")
    print(f"\n  Rankings saved → {tsv_path}")


def main():
    ap = argparse.ArgumentParser(description='dMaSIF-inspired C-B interaction screen')
    ap.add_argument('--candidates', nargs='*',
                    help='PDB paths to screen as protein C (default: new downloads)')
    ap.add_argument('--clear-cache', action='store_true',
                    help='Force recomputation of A-interface fingerprints')
    ap.add_argument('--use-esmfold', action='store_true',
                    help='Build reference fingerprint from ESMFold standalone A '
                         f'(pLDDT>{PLDDT_THRESH:.0f} filtered) instead of bound complex')
    ap.add_argument('--use-apo', action='store_true',
                    help=f'Build reference fingerprint from apo crystal monomer ({APO_A_PDB}), '
                         f'excluding induced-fit residues {sorted(EXCLUDE_RNUMS)}')
    ap.add_argument('--use-sgk', action='store_true',
                    help=f'Build reference fingerprint from {SGK_A_PDB} (apo DT), '
                         f'NAD-site residues identified via superposition onto {F0L_A_PDB}')
    ap.add_argument('--rankings-tsv', default='results/dmasif/rankings.tsv',
                    help='Output TSV path for rankings (default: results/dmasif/rankings.tsv)')
    args = ap.parse_args()

    cache = A_CACHE_SGK if args.use_sgk else (A_CACHE_APO if args.use_apo else (A_CACHE_ESM if args.use_esmfold else A_CACHE))
    if args.clear_cache and os.path.exists(cache):
        os.remove(cache)
        print(f"Removed cache: {cache}")

    t0         = time.time()
    pdb_parser = PDBParser(QUIET=True)

    A_fps, A_iface_coords, B_pos, B_vdw = load_or_build_A_fingerprints(
        pdb_parser, use_esmfold=args.use_esmfold, use_apo=args.use_apo, use_sgk=args.use_sgk)
    print(f"  A interface fingerprints : {len(A_fps)}")
    print(f"  Chain B surface atoms    : {len(B_pos)}")

    if args.candidates:
        candidates = args.candidates
    else:
        import glob
        candidates = sorted(glob.glob('binding/*.pdb'))

    results = screen_candidates(candidates, A_fps, A_iface_coords, B_pos, B_vdw, pdb_parser)
    print_summary(results, tsv_path=args.rankings_tsv)

    print(f"Total runtime: {time.time()-t0:.0f}s")
    print(f"Docked PDB files saved in: results/dmasif/")


if __name__ == '__main__':
    main()
