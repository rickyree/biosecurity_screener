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

import os, sys, copy, time, warnings, argparse, subprocess
import numpy as np
from scipy.spatial import KDTree
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.SASA import ShrakeRupley

AMINA = os.path.join(os.path.dirname(sys.executable), 'amina')

warnings.filterwarnings('ignore')
np.random.seed(42)

os.makedirs('results/dmasif', exist_ok=True)

SGK_A_PDB            = 'binding/1sgk_A.pdb'   # default apo receptor
DEFAULT_BOUND_REF    = 'binding/1f0l_A.pdb'   # default holo reference (chain A + APU)
DEFAULT_BOUND_LIGAND = 'APU'                   # default ligand residue name

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


def resolve_bound_chain(structure, bound_chain):
    """Return (receptor_chain_id, partner_chain_id_or_None).

    1 protein chain  → receptor=that chain, partner=None (use ligand for interface)
    2 protein chains → receptor=chain A if present else first, partner=the other chain
    >2 protein chains → error (ambiguous multi-partner complex)

    If bound_chain is explicitly provided it is used as the receptor chain.
    """
    protein_chains = [
        chain.id
        for model in structure
        for chain in model
        if any('CA' in res for res in chain if res.id[0] == ' ')
    ]

    if len(protein_chains) > 2:
        raise SystemExit(
            f"ERROR: bound reference contains {len(protein_chains)} protein chains "
            f"({protein_chains}). Multi-partner complexes are not supported — the "
            f"workflow cannot distinguish between distinct interfaces. Provide a "
            f"single-chain or two-chain (receptor + one partner) structure as --bound-ref."
        )

    if bound_chain:
        all_chains = [chain.id for model in structure for chain in model]
        if bound_chain not in all_chains:
            raise SystemExit(
                f"ERROR: --bound-chain '{bound_chain}' not found in bound reference "
                f"(chains present: {all_chains})."
            )
        receptor = bound_chain
    elif len(protein_chains) == 0:
        all_chains = [chain.id for model in structure for chain in model]
        if 'A' in all_chains:
            receptor = 'A'
        else:
            raise SystemExit(
                f"ERROR: no protein chains detected in bound reference "
                f"(chains: {all_chains}). Use --bound-chain to specify."
            )
    elif len(protein_chains) == 1:
        receptor = protein_chains[0]
    else:  # 2 protein chains
        receptor = 'A' if 'A' in protein_chains else protein_chains[0]

    partner = None
    if len(protein_chains) == 2:
        partner = next(c for c in protein_chains if c != receptor)

    return receptor, partner


def get_bound_partner_coords(structure, receptor_chain, partner_chain, bound_ligand):
    """Return coords of the interface-defining atoms from the bound reference.

    Two-chain complex: return all heavy atoms of partner_chain.
    Single-chain:      return atoms of the ligand residue (bound_ligand name).
    """
    if partner_chain is not None:
        coords = [a.coord for model in structure for chain in model
                  if chain.id == partner_chain
                  for res in chain if res.id[0] == ' '
                  for a in res
                  if (a.element or a.name[0]).strip().upper() != 'H']
        if not coords:
            raise SystemExit(
                f"ERROR: no heavy atoms found in partner chain '{partner_chain}' "
                f"of bound reference."
            )
        return np.array(coords)
    else:
        if not bound_ligand:
            raise SystemExit(
                "ERROR: --bound-ref has a single protein chain so a ligand is needed "
                "to define the interface. Specify it with --bound-ligand <RES>."
            )
        coords = [a.coord for model in structure for chain in model
                  if chain.id == receptor_chain
                  for res in chain
                  if res.get_resname().strip() == bound_ligand
                  for a in res]
        if not coords:
            raise SystemExit(
                f"ERROR: ligand '{bound_ligand}' not found in chain '{receptor_chain}' "
                f"of bound reference. Use --bound-ligand to specify the correct residue name."
            )
        return np.array(coords)


def get_union_iface_res(apo_ref, bound_ref=None, bound_chain=None, bound_ligand=None,
                        pesto_threshold=0.5, work_dir='results/dmasif/union_tmp',
                        pesto_type=None):
    """Return sorted list of apo residue numbers from PeSTo ∪ P2Rank pocket-1.

    PeSTo interface type is chosen automatically when bound_ref is provided:
      - 1 protein chain in bound_ref → Protein-Ligand-Interface
      - 2 protein chains             → Protein-Protein-Interface
    If bound_ref is None, pesto_type must be 'ligand' or 'protein' (default: 'ligand').
    """
    import pandas as pd

    os.makedirs(work_dir, exist_ok=True)
    if bound_ref is not None:
        pdb_parser = PDBParser(QUIET=True)
        bound_struct = pdb_parser.get_structure('BOUND', bound_ref)
        _, partner_chain = resolve_bound_chain(bound_struct, bound_chain)
        iface_type = ('Protein-Protein-Interface' if partner_chain
                      else 'Protein-Ligand-Interface')
        print(f"  [union] Partner type: {'protein chain ' + partner_chain if partner_chain else 'ligand ' + str(bound_ligand)}")
    else:
        explicit = pesto_type or 'ligand'
        iface_type = ('Protein-Protein-Interface' if explicit == 'protein'
                      else 'Protein-Ligand-Interface')
        print(f"  [union] No bound reference — using PeSTo type: {iface_type}")
    print(f"  [union] PeSTo interface type: {iface_type}")

    # Step 1: clean apo with pdb-cleaner
    apo_stem   = os.path.splitext(os.path.basename(apo_ref))[0]
    clean_dir  = os.path.join(work_dir, 'cleaned')
    clean_path = os.path.join(clean_dir, f'{apo_stem}_cleaned.pdb')
    os.makedirs(clean_dir, exist_ok=True)
    if not os.path.exists(clean_path):
        print(f"  [union] Cleaning apo with pdb-cleaner …", flush=True)
        subprocess.run([AMINA, 'run', 'pdb-cleaner', '--pdb', apo_ref,
                        '-o', clean_dir, '-j', apo_stem],
                       capture_output=True)
    if not os.path.exists(clean_path):
        print(f"  [union] pdb-cleaner failed; using original apo.", flush=True)
        clean_path = apo_ref
    else:
        print(f"  [union] Using cleaned apo: {clean_path}")

    # Step 2: run PeSTo and P2Rank (skip if results already cached)
    pesto_dir  = os.path.join(work_dir, 'pesto')
    p2rank_dir = os.path.join(work_dir, 'p2rank')
    os.makedirs(pesto_dir,  exist_ok=True)
    os.makedirs(p2rank_dir, exist_ok=True)

    pesto_csv_path  = os.path.join(pesto_dir,
                                   f'pesto_{apo_stem}_bfactor_{iface_type}_residues.csv')
    p2rank_csv_path = os.path.join(p2rank_dir, f'p2rank_{apo_stem}_residues.csv')

    def submit(args):
        r = subprocess.run([AMINA, 'run'] + args + ['--background'],
                           capture_output=True, text=True)
        for line in (r.stdout + r.stderr).splitlines():
            if 'Job submitted:' in line:
                return line.split('Job submitted:')[1].strip()
        return None

    jobs = []
    if os.path.exists(pesto_csv_path):
        print(f"  [union] PeSTo results already cached, skipping.", flush=True)
    else:
        print(f"  [union] Submitting PeSTo ({iface_type}) …", flush=True)
        job = submit(['pesto', '--pdb', clean_path,
                      '--threshold', str(pesto_threshold),
                      '-o', pesto_dir, '-j', apo_stem])
        if job:
            jobs.append((job, 'PeSTo', pesto_dir))

    if os.path.exists(p2rank_csv_path):
        print(f"  [union] P2Rank results already cached, skipping.", flush=True)
    else:
        print(f"  [union] Submitting P2Rank …", flush=True)
        job = submit(['p2rank', '--pdb', clean_path,
                      '-o', p2rank_dir, '-j', apo_stem])
        if job:
            jobs.append((job, 'P2Rank', p2rank_dir))

    # Step 3: wait and download only new jobs
    for job_id, label, out_dir in jobs:
        print(f"  [union] Waiting for {label} (job {job_id[:8]}) …", flush=True)
        subprocess.run([AMINA, 'jobs', 'wait', job_id], check=True)
        subprocess.run([AMINA, 'jobs', 'download', job_id, '-o', out_dir], check=True)

    # Step 4: parse PeSTo residues for the chosen interface type
    pesto_res = set()
    if os.path.exists(pesto_csv_path):
        df = pd.read_csv(pesto_csv_path)
        pesto_res = set(df['residue_number'].tolist())
        print(f"  [union] PeSTo {iface_type}: {len(pesto_res)} residues")
    else:
        print(f"  [union] WARNING: PeSTo CSV not found at {pesto_csv_path}", file=sys.stderr)

    # Step 5: parse P2Rank pocket-1 residues only
    p2rank_res = set()
    if os.path.exists(p2rank_csv_path):
        df = pd.read_csv(p2rank_csv_path, skipinitialspace=True)
        p2rank_res = set(df[df['pocket'] == 1]['residue_label'].tolist())
        print(f"  [union] P2Rank pocket-1: {len(p2rank_res)} residues")
    else:
        print(f"  [union] WARNING: P2Rank CSV not found at {p2rank_csv_path}", file=sys.stderr)

    union = sorted(pesto_res | p2rank_res)
    print(f"  [union] Union: {len(union)} residues → {union}")
    return union


def load_or_build_A_fingerprints(pdb_parser, apo_ref, bound_ref=None,
                                  bound_chain=None, bound_ligand=None,
                                  plddt_thresh=0.0, use_union=False,
                                  pesto_type=None, iface_residues=None):
    """Return (A_fps, A_iface_coords, B_pos, B_vdw) either from cache or by computing.

    apo_ref:        path to the apo receptor PDB.
    bound_ref:      path to holo PDB defining the binding site. Required unless
                    use_union=True or iface_residues is provided.
    bound_chain:    chain in bound_ref containing the receptor (auto-detected).
    bound_ligand:   residue name of the ligand in bound_ref (single-chain case).
    plddt_thresh:   drop interface residues whose CA B-factor < this value.
    use_union:      derive interface via PeSTo+P2Rank union instead of <5Å.
    pesto_type:     'ligand' or 'protein' — used when use_union=True and
                    bound_ref is None.
    iface_residues: explicit list of residue numbers to use as the interface,
                    bypassing both <5Å and union methods. bound_ref not needed.
    """
    apo_stem   = os.path.splitext(os.path.basename(apo_ref))[0]
    bound_stem = os.path.splitext(os.path.basename(bound_ref))[0] if bound_ref else 'noref'
    if iface_residues is not None:
        suffix = '_explicit'
    elif use_union:
        suffix = '_union'
    else:
        suffix = ''
    cache = f'results/dmasif_cache/A_interface_cache_{apo_stem}_{bound_stem}{suffix}.npz'
    if os.path.exists(cache):
        print(f"  Loading cached A-interface fingerprints from {cache}")
        d = np.load(cache)
        return d['A_fps'], d['A_iface_coords'], d['B_pos'], d['B_vdw']

    if iface_residues is not None:
        method_str = f"explicit residues {sorted(iface_residues)}"
    elif use_union:
        method_str = "PeSTo+P2Rank union"
    else:
        if not bound_ref or not os.path.exists(bound_ref):
            sys.exit(f"ERROR: --bound-ref is required when not using --union or --iface-residues")
        bound_struct = pdb_parser.get_structure('BOUND', bound_ref)
        receptor_chain, partner_chain = resolve_bound_chain(bound_struct, bound_chain)
        iface_label = (f"chain {partner_chain} contact" if partner_chain
                       else f"ligand {bound_ligand}")
        method_str = f"<5Å from {iface_label}"

    plddt_str  = f", pLDDT>{plddt_thresh:.0f}" if plddt_thresh > 0 else ""
    print("=" * 65)
    print(f" Step 1: Building dMaSIF fingerprint for A's interface")
    print(f"         ({os.path.basename(apo_ref)}, {method_str}{plddt_str})")
    print("=" * 65)

    A_fps, A_fp_coords = [], []
    B_pos = np.empty((0, 3))
    B_vdw = np.empty(0)

    if not os.path.exists(apo_ref):
        sys.exit(f"ERROR: apo reference not found at {apo_ref}")
    apo_struct   = pdb_parser.get_structure('APO', apo_ref)
    apo_residues = [res for model in apo_struct for chain in model
                    for res in chain if res.id[0] == ' ' and 'CA' in res]
    apo_label = os.path.basename(apo_ref)

    if iface_residues is not None:
        # ── Explicit residues path ────────────────────────────────────────────
        iface_res = [res for res in apo_residues if res.id[1] in set(iface_residues)]
        print(f"  Interface residues explicit ({len(iface_res)}): "
              f"{[r.id[1] for r in iface_res]}")
    elif use_union:
        # ── Union path: PeSTo + P2Rank ────────────────────────────────────────
        union_resnums = set(get_union_iface_res(
            apo_ref, bound_ref, bound_chain, bound_ligand))
        iface_res = [res for res in apo_residues if res.id[1] in union_resnums]
        print(f"  Interface residues from union ({len(iface_res)}): "
              f"{[r.id[1] for r in iface_res]}")
    else:
        # ── 5Å path: superpose apo onto holo, distance cutoff ─────────────────
        from Bio.SeqUtils import seq1
        from Bio import pairwise2

        bound_residues = [res for model in bound_struct for chain in model
                          if chain.id == receptor_chain
                          for res in chain if res.id[0] == ' ' and 'CA' in res]
        bound_ca  = np.array([res['CA'].coord for res in bound_residues])
        bound_seq = ''.join(seq1(r.get_resname()) for r in bound_residues)

        partner_coords = get_bound_partner_coords(
            bound_struct, receptor_chain, partner_chain, bound_ligand)
        print(f"  Interface-defining atoms in bound ref: {len(partner_coords)}")

        apo_ca  = np.array([res['CA'].coord for res in apo_residues])
        apo_seq = ''.join(seq1(r.get_resname()) for r in apo_residues)

        aln = pairwise2.align.globalds(
            apo_seq, bound_seq,
            pairwise2.substitution_matrices.load("BLOSUM62"),
            -10, -0.5, one_alignment_only=True)[0]
        apo_aln, bnd_aln = aln.seqA, aln.seqB
        apo_idx, bnd_idx = [], []
        si, fi = 0, 0
        for a, b in zip(apo_aln, bnd_aln):
            if a != '-' and b != '-':
                apo_idx.append(si)
                bnd_idx.append(fi)
            if a != '-': si += 1
            if b != '-': fi += 1
        apo_ca_aln   = apo_ca[apo_idx]
        bound_ca_aln = bound_ca[bnd_idx]

        R_sup, t_sup = superimpose_kabsch(apo_ca_aln, bound_ca_aln)
        rmsd_sup = float(np.sqrt(np.mean(
            np.sum(((R_sup @ apo_ca_aln.T).T + t_sup - bound_ca_aln)**2, axis=1))))
        print(f"  Superimposed {apo_label} onto {os.path.basename(bound_ref)} "
              f"chain {receptor_chain} ({len(apo_idx)} Cα pairs, RMSD={rmsd_sup:.2f} Å)")

        # Identify interface residues directly in holo frame, then map to apo
        # via sequence alignment — avoids querying apo atoms in a transformed frame.
        partner_tree = KDTree(partner_coords)
        holo_iface_idx = []
        for i, res in enumerate(bound_residues):
            for atom in res:
                if partner_tree.query_ball_point(atom.coord, r=5.0):
                    holo_iface_idx.append(i)
                    break

        bnd_to_apo = {bi: ai for ai, bi in zip(apo_idx, bnd_idx)}
        iface_res = []
        for hi in holo_iface_idx:
            if hi in bnd_to_apo:
                iface_res.append(apo_residues[bnd_to_apo[hi]])

    if plddt_thresh > 0:
        iface_res = [res for res in iface_res
                     if res['CA'].get_bfactor() >= plddt_thresh]
        print(f"  Interface residues after pLDDT>{plddt_thresh:.0f} filter: {len(iface_res)}")
    elif not use_union:
        print(f"  Interface residues ({len(iface_res)}): {[r.id[1] for r in iface_res]}")

    A_atoms, A_pos = all_heavy_atoms(apo_struct)
    print(f"  {apo_label} heavy atoms: {len(A_atoms)}")
    A_desc = build_desc(A_atoms, A_pos)
    A_tree = KDTree(A_pos)

    for res in iface_res:
        if 'CA' not in res:
            continue
        ca_coord = res['CA'].coord
        _, idx = A_tree.query(ca_coord)
        fp = patch_fp(idx, A_pos, A_desc, A_tree)
        if fp is not None:
            A_fps.append(fp)
            A_fp_coords.append(A_pos[idx])

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

        out_pdb = f"results/dmasif/structures/{name}_docked.pdb"
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
    ap.add_argument('--apo-ref', required=True, metavar='PDB',
                    help='Apo receptor PDB to build interface fingerprints from.')
    ap.add_argument('--plddt-thresh', type=float, default=0.0, metavar='THRESH',
                    help='Drop interface residues with CA B-factor below this value '
                         '(default: 0 = disabled). Set >0 when --apo-ref is a '
                         'predicted structure storing pLDDT in the B-factor column.')
    ap.add_argument('--bound-ref', default=None, metavar='PDB',
                    help='Holo PDB defining the binding site. '
                         '1 protein chain: ligand proximity defines interface. '
                         '2 protein chains: chain-chain interface used directly. '
                         'Not required when --union or --iface-residues is used.')
    ap.add_argument('--bound-chain', default=None, metavar='ID',
                    help='Chain ID of the receptor in --bound-ref (auto-detected if omitted).')
    ap.add_argument('--bound-ligand', default=None, metavar='RES',
                    help='Residue name of the ligand in --bound-ref (single-chain case only). '
                         'Required when --bound-ref has one protein chain.')
    ap.add_argument('--rankings-tsv', default=None,
                    help='Output TSV path for rankings. Defaults to '
                         'results/dmasif/rankings_union.tsv (--union) or '
                         'results/dmasif/rankings.tsv (default).')
    ap.add_argument('--union', action='store_true',
                    help='Derive interface residues via PeSTo+P2Rank union '
                         'instead of <5Å from holo. Uses a separate cache and '
                         'outputs rankings to rankings_union.tsv by default.')
    ap.add_argument('--iface-residues', default=None, metavar='NUMS',
                    help='Comma-separated residue numbers to use as the binding site '
                         '(e.g. 32,64,65,68,125,221). Bypasses <5Å and union methods; '
                         '--bound-ref not required.')
    args = ap.parse_args()

    iface_residues = ([int(x) for x in args.iface_residues.split(',')]
                      if args.iface_residues else None)

    t0         = time.time()
    pdb_parser = PDBParser(QUIET=True)

    # Derive cache path so different apo / bound-ref / method combos don't collide
    apo_stem   = os.path.splitext(os.path.basename(args.apo_ref))[0]
    bound_stem = os.path.splitext(os.path.basename(args.bound_ref))[0] if args.bound_ref else 'noref'
    if iface_residues is not None:
        suffix = '_explicit'
    elif args.union:
        suffix = '_union'
    else:
        suffix = ''
    cache  = f'results/dmasif_cache/A_interface_cache_{apo_stem}_{bound_stem}{suffix}.npz'
    if args.clear_cache and os.path.exists(cache):
        os.remove(cache)
        print(f"Removed cache: {cache}")

    if args.rankings_tsv is None:
        if iface_residues is not None:
            args.rankings_tsv = f'results/dmasif/rankings_{apo_stem}_explicit.tsv'
        elif args.union:
            args.rankings_tsv = f'results/dmasif/rankings_{apo_stem}_union.tsv'
        else:
            args.rankings_tsv = 'results/dmasif/rankings.tsv'

    A_fps, A_iface_coords, B_pos, B_vdw = load_or_build_A_fingerprints(
        pdb_parser, apo_ref=args.apo_ref,
        bound_ref=args.bound_ref,
        bound_chain=args.bound_chain,
        bound_ligand=args.bound_ligand,
        plddt_thresh=args.plddt_thresh,
        use_union=args.union,
        iface_residues=iface_residues)
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
    print(f"Docked PDB files saved in: results/dmasif/structures/")


if __name__ == '__main__':
    main()
