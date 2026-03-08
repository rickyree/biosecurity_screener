#!/usr/bin/env python3
"""
Standalone ESM2 pre-filter for candidate sequences with integrated interface identification.

Two filters (both position-free: reference binding site embeddings are matched
to candidate positions via nearest-neighbour search in embedding space):
  1. Inter-cluster contact score: ESM2 contact probability between candidate
     positions that best match the reference binding site cluster embeddings.
  2. Chemistry similarity: mean cosine similarity of reference binding site
     embeddings to their nearest neighbours in the candidate patch.

Binding site identification integrated directly (no dMaSIF cache dependency):
  - <5Å method: residues within 5Å of ligand/partner (default for holo structures)
  - PeSTo method: ML-based binding site predictions (automatic for --apo structures)
  - Explicit method: user-provided residue numbers (overrides automatic selection)
"""
#final file
import sys, os, argparse, hashlib, subprocess
import numpy as np
import pandas as pd
import torch
import esm as esmlib
from scipy.spatial import KDTree
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import pairwise2

AMINA = os.path.join(os.path.dirname(sys.executable), 'amina')
ESM2_CACHE_DIR = 'results/esm2_embed_cache'

# ── ESM2 inference (with sequence-keyed cache) ────────────────────────────────

def seq_hash(seq):
    """SHA256 hash of the cleaned sequence — used as cache key."""
    clean = seq.replace('X', 'G').replace('-', 'G')[:1000]
    return hashlib.sha256(clean.encode()).hexdigest()

def get_outputs(model, batch_converter, seq, label):
    clean = seq.replace('X', 'G').replace('-', 'G')[:1000]
    cache_path = os.path.join(ESM2_CACHE_DIR, f'{seq_hash(seq)}.npz')

    if os.path.exists(cache_path):
        d = np.load(cache_path)
        return d['emb'], d['contacts']

    data = [(label, clean)]
    _, _, tokens = batch_converter(data)
    device = next(model.parameters()).device
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[33], return_contacts=True)
    emb = out['representations'][33][0, 1:1+len(clean)].cpu().numpy()
    contacts = out['contacts'][0].cpu().numpy()

    os.makedirs(ESM2_CACHE_DIR, exist_ok=True)
    np.savez(cache_path, emb=emb, contacts=contacts)
    return emb, contacts

# ── Interface identification functions ────────────────────────────────────────

def resolve_bound_chain(structure, bound_chain):
    """Return (receptor_chain_id, partner_chain_id_or_None).

    1 protein chain  → receptor=that chain, partner=None (use ligand for interface)
    2 protein chains → receptor=chain A if present else first, partner=the other chain
    >2 protein chains → error (ambiguous multi-partner complex)
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
            f"({protein_chains}). Multi-partner complexes are not supported."
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

def superimpose_kabsch(mobile_ca, target_ca):
    """Return R, t that maps mobile CA coords onto target CA coords (Kabsch)."""
    cp, cq = mobile_ca.mean(0), target_ca.mean(0)
    Pc, Qc = mobile_ca - cp, target_ca - cq
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp

def identify_interface_5a(apo_ref, bound_ref, bound_chain=None, bound_ligand=None, pdb_parser=None):
    """Identify interface residues using <5Å distance method."""
    if pdb_parser is None:
        pdb_parser = PDBParser(QUIET=True)

    bound_struct = pdb_parser.get_structure('BOUND', bound_ref)
    receptor_chain, partner_chain = resolve_bound_chain(bound_struct, bound_chain)

    partner_coords = get_bound_partner_coords(
        bound_struct, receptor_chain, partner_chain, bound_ligand)
    print(f"  Interface-defining atoms in bound ref: {len(partner_coords)}")

    apo_struct = pdb_parser.get_structure('APO', apo_ref)
    apo_residues = [res for model in apo_struct for chain in model
                    for res in chain if res.id[0] == ' ' and 'CA' in res]

    bound_residues = [res for model in bound_struct for chain in model
                      if chain.id == receptor_chain
                      for res in chain if res.id[0] == ' ' and 'CA' in res]

    # Sequence alignment and superposition
    bound_ca = np.array([res['CA'].coord for res in bound_residues])
    bound_seq = ''.join(seq1(r.get_resname()) for r in bound_residues)

    apo_ca = np.array([res['CA'].coord for res in apo_residues])
    apo_seq = ''.join(seq1(r.get_resname()) for r in apo_residues)

    aln = pairwise2.align.globalds(
        apo_seq, bound_seq,
        pairwise2.substitution_matrices.load("BLOSUM62"),
        -10, -0.5, one_alignment_only=True)[0]

    apo_idx, bnd_idx = [], []
    si, fi = 0, 0
    for a, b in zip(aln.seqA, aln.seqB):
        if a != '-' and b != '-':
            apo_idx.append(si)
            bnd_idx.append(fi)
        if a != '-': si += 1
        if b != '-': fi += 1

    apo_ca_aln = apo_ca[apo_idx]
    bound_ca_aln = bound_ca[bnd_idx]

    R_sup, t_sup = superimpose_kabsch(apo_ca_aln, bound_ca_aln)
    rmsd_sup = float(np.sqrt(np.mean(
        np.sum(((R_sup @ apo_ca_aln.T).T + t_sup - bound_ca_aln)**2, axis=1))))
    print(f"  Superimposed apo onto bound chain {receptor_chain} "
          f"({len(apo_idx)} Cα pairs, RMSD={rmsd_sup:.2f} Å)")

    # Identify interface residues
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

    return iface_res

def run_pesto_prediction(apo_ref, bound_ref, bound_chain=None, bound_ligand=None,
                        pesto_threshold=0.5, work_dir='results/interface_tmp'):
    """Run PeSTo prediction method.

    Uses only the highest-confidence pocket for improved specificity.
    Falls back to threshold-based filtering if pocket information unavailable.
    """
    os.makedirs(work_dir, exist_ok=True)
    pdb_parser = PDBParser(QUIET=True)

    bound_struct = pdb_parser.get_structure('BOUND', bound_ref) if bound_ref else None
    if bound_struct:
        _, partner_chain = resolve_bound_chain(bound_struct, bound_chain)
        iface_type = ('Protein-Protein-Interface' if partner_chain
                      else 'Protein-Ligand-Interface')
    else:
        iface_type = 'Protein-Ligand-Interface'  # default

    print(f"  [pesto] PeSTo interface type: {iface_type}")

    # Clean apo with pdb-cleaner
    apo_stem = os.path.splitext(os.path.basename(apo_ref))[0]
    clean_dir = os.path.join(work_dir, 'cleaned')
    clean_path = os.path.join(clean_dir, f'{apo_stem}_cleaned.pdb')
    os.makedirs(clean_dir, exist_ok=True)

    if not os.path.exists(clean_path):
        print(f"  [pesto] Cleaning apo with pdb-cleaner …", flush=True)
        subprocess.run([AMINA, 'run', 'pdb-cleaner', '--pdb', apo_ref,
                        '-o', clean_dir, '-j', apo_stem],
                       capture_output=True)
    if not os.path.exists(clean_path):
        print(f"  [pesto] pdb-cleaner failed; using original apo.", flush=True)
        clean_path = apo_ref

    # Submit job
    def submit(args):
        r = subprocess.run([AMINA, 'run'] + args + ['--background'],
                           capture_output=True, text=True)
        for line in (r.stdout + r.stderr).splitlines():
            if 'Job submitted:' in line:
                return line.split('Job submitted:')[1].strip()
        return None

    pesto_dir = os.path.join(work_dir, 'pesto')
    os.makedirs(pesto_dir, exist_ok=True)

    pesto_csv_path = os.path.join(pesto_dir,
                                  f'pesto_{apo_stem}_bfactor_{iface_type}_residues.csv')

    job_id = None
    if not os.path.exists(pesto_csv_path):
        print(f"  [pesto] Submitting PeSTo ({iface_type}) …", flush=True)
        job_id = submit(['pesto', '--pdb', clean_path,
                         '--threshold', str(pesto_threshold),
                         '-o', pesto_dir, '-j', apo_stem])

    # Wait and download
    if job_id:
        print(f"  [pesto] Waiting for PeSTo (job {job_id[:8]}) …", flush=True)
        subprocess.run([AMINA, 'jobs', 'wait', job_id], check=True)
        subprocess.run([AMINA, 'jobs', 'download', job_id, '-o', pesto_dir], check=True)

    # Parse results
    pesto_res = set()
    if os.path.exists(pesto_csv_path):
        df = pd.read_csv(pesto_csv_path)

        # Use only the highest confidence pocket for specificity
        if 'pocket' in df.columns or 'pocket_id' in df.columns:
            pocket_col = 'pocket_id' if 'pocket_id' in df.columns else 'pocket'

            # Find pocket with highest mean binding probability
            if 'binding_probability' in df.columns:
                pocket_scores = df.groupby(pocket_col)['binding_probability'].mean()
                top_pocket = pocket_scores.idxmax()
                print(f"  [pesto] Using top pocket {top_pocket} (score: {pocket_scores[top_pocket]:.3f})")
            else:
                # Fallback: use pocket 1 if no binding probability column
                top_pocket = 1
                print(f"  [pesto] Using pocket {top_pocket} (no scores available)")

            # Filter to top pocket only
            df_filtered = df[df[pocket_col] == top_pocket]
            pesto_res = set(df_filtered['residue_number'].tolist())
            print(f"  [pesto] Top pocket {top_pocket}: {len(pesto_res)} residues")
        else:
            # Fallback: use all residues above threshold if no pocket information
            if 'binding_probability' in df.columns:
                df_filtered = df[df['binding_probability'] >= pesto_threshold]
                pesto_res = set(df_filtered['residue_number'].tolist())
                print(f"  [pesto] All pockets (threshold ≥{pesto_threshold}): {len(pesto_res)} residues")
            else:
                pesto_res = set(df['residue_number'].tolist())
                print(f"  [pesto] All residues: {len(pesto_res)} residues")

    pesto_residues = sorted(pesto_res)
    print(f"  [pesto] Selected residues: {pesto_residues}")
    return pesto_residues

def derive_binding_indices_direct(structure_ref, bound_ref=None, bound_chain=None,
                                 bound_ligand=None, use_pesto=False,
                                 iface_residues=None, plddt_thresh=0.0, gap=10):
    """Derive binding site indices directly without requiring dMaSIF cache.

    Interface identification methods:
    - iface_residues: explicit list of residue numbers
    - use_pesto: PeSTo predictions on apo structure
    - default: <5Å from bound partner/ligand in holo structure
    """
    pdb_parser = PDBParser(QUIET=True)
    struct = pdb_parser.get_structure('STRUCT', structure_ref)
    struct_residues = [res for model in struct for chain in model
                       for res in chain if res.id[0] == ' ' and 'CA' in res]

    if iface_residues is not None:
        print(f"  Using explicit interface residues: {sorted(iface_residues)}")
        iface_res = [res for res in struct_residues if res.id[1] in set(iface_residues)]
    elif use_pesto:
        print("  Using PeSTo method")
        pesto_resnums = set(run_pesto_prediction(
            structure_ref, bound_ref, bound_chain, bound_ligand))
        iface_res = [res for res in struct_residues if res.id[1] in pesto_resnums]
    else:
        print(f"  Using <5Å method from {os.path.basename(structure_ref)}")
        iface_res = identify_interface_5a(structure_ref, structure_ref, bound_chain, bound_ligand, pdb_parser)

    # pLDDT filtering
    if plddt_thresh > 0:
        iface_res = [res for res in iface_res
                     if res['CA'].get_bfactor() >= plddt_thresh]
        print(f"  Interface residues after pLDDT>{plddt_thresh:.0f} filter: {len(iface_res)}")

    # Convert to sequence indices
    res_to_idx = {res.id[1]: i for i, res in enumerate(struct_residues)}
    seq_idx = sorted([res_to_idx[res.id[1]] for res in iface_res if res.id[1] in res_to_idx])

    # Cluster by sequence gaps
    if not seq_idx:
        sys.exit("ERROR: no valid interface residues found")

    clusters, current = [], [seq_idx[0]]
    for i in range(1, len(seq_idx)):
        if seq_idx[i] - seq_idx[i-1] >= gap:
            clusters.append(current)
            current = []
        current.append(seq_idx[i])
    clusters.append(current)

    print(f"  Binding site: {len(seq_idx)} residues, "
          f"{len(clusters)} cluster(s) starting at: {[c[0] for c in clusters]}")
    print(f"  Interface residues: {[struct_residues[i].id[1] for i in seq_idx]}")

    return seq_idx, clusters

# ── Matching and scoring ──────────────────────────────────────────────────────

def match_clusters(ref_emb, cand_emb, cand_contacts, ref_cluster_idx, contact_thresh=0.1):
    """For each reference cluster index, find the nearest-neighbour position
    in the candidate via cosine similarity in embedding space."""
    cand_norm = cand_emb / (np.linalg.norm(cand_emb, axis=1, keepdims=True) + 1e-8)

    matched = []
    for ri in ref_cluster_idx:
        if ri >= len(ref_emb):
            continue
        rv = ref_emb[ri]
        rv_norm = rv / (np.linalg.norm(rv) + 1e-8)
        matched.append(int(np.argmax(cand_norm @ rv_norm)))

    if not matched:
        return []

    matched = list(dict.fromkeys(matched))   # deduplicate

    if len(matched) == 1:
        return matched

    # Contact-based compaction
    L = cand_contacts.shape[0]
    valid = [i for i in matched if i < L]
    if len(valid) < 2:
        return valid

    sub = cand_contacts[np.ix_(valid, valid)]
    np.fill_diagonal(sub, 0.0)
    mean_contact = sub.mean(axis=1)
    compact = [valid[i] for i in range(len(valid))
               if mean_contact[i] >= contact_thresh]

    return compact if compact else valid

def contact_score_matched(contacts, cluster_matches):
    """Inter-cluster contact score across all pairs of matched cluster sets."""
    L = contacts.shape[0]
    valid = [[i for i in m if i < L] for m in cluster_matches]
    pairs = []
    for a in range(len(valid)):
        for b in range(a + 1, len(valid)):
            if valid[a] and valid[b]:
                pairs.append(contacts[np.ix_(valid[a], valid[b])].mean())
    return float(np.mean(pairs)) if pairs else float('nan')

def chem_sim_matched(ref_emb, cand_emb, ref_binding_idx, matched_cand_idx):
    """Mean cosine similarity between reference binding site embeddings and
    their best match within the already-matched, compacted candidate patch."""
    if len(matched_cand_idx) < 3:
        return float('nan')

    patch = cand_emb[matched_cand_idx]
    patch_norm = patch / (np.linalg.norm(patch, axis=1, keepdims=True) + 1e-8)

    sims = []
    for ri in ref_binding_idx:
        if ri >= len(ref_emb):
            continue
        rv = ref_emb[ri]
        rv_norm = rv / (np.linalg.norm(rv) + 1e-8)
        sims.append(float((patch_norm @ rv_norm).max()))

    if len(sims) < 3:
        return float('nan')
    return float(np.mean(sims))

def ref_contact_score(ref_contacts, clusters):
    """Contact score on the reference itself using its own cluster indices."""
    L = ref_contacts.shape[0]
    valid = [[i for i in c if i < L] for c in clusters]
    pairs = []
    for a in range(len(valid)):
        for b in range(a + 1, len(valid)):
            if valid[a] and valid[b]:
                pairs.append(ref_contacts[np.ix_(valid[a], valid[b])].mean())
    return float(np.mean(pairs)) if pairs else float('nan')

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--structure', required=True, metavar='PDB',
                    help='Input protein structure (apo or holo).')
    ap.add_argument('--apo', action='store_true',
                    help='Structure is apo (unbound) - automatically uses PeSTo for binding site prediction.')
    ap.add_argument('--bound-chain', default=None, metavar='ID',
                    help='Chain ID of the receptor in holo structure (auto-detected if omitted).')
    ap.add_argument('--bound-ligand', default=None, metavar='RES',
                    help='Residue name of the ligand in holo structure (single-chain case only).')
    ap.add_argument('--candidates-tsv', required=True, metavar='TSV',
                    help='TSV file of candidate sequences with columns "candidate" and "sequence_chain1".')
    ap.add_argument('--out', default=None, metavar='TSV',
                    help='Output TSV path (default: output_{structure_stem}[_apo|_explicit].tsv)')
    ap.add_argument('--iface-residues', default=None, metavar='NUMS',
                    help='Comma-separated residue numbers to use as the binding site (e.g. 32,64,65).')
    ap.add_argument('--plddt-thresh', type=float, default=0.0, metavar='THRESH',
                    help='Drop interface residues with CA B-factor below this value (default: 0).')
    args = ap.parse_args()

    iface_residues_list = ([int(x) for x in args.iface_residues.split(',')]
                          if args.iface_residues else None)

    structure_stem = os.path.splitext(os.path.basename(args.structure))[0]
    if iface_residues_list is not None:
        suffix = '_explicit'
    elif args.apo:
        suffix = '_apo'
    else:
        suffix = ''

    out_path = args.out or f'output_{structure_stem}{suffix}.tsv'

    print("Loading ESM2 650M...", flush=True)
    model, alphabet = esmlib.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model = model.to(device)
    model.eval()
    print(f"  device: {device}", flush=True)

    # ── Interface identification ──
    print(f"Interface identification for {args.structure}", flush=True)

    # Determine method based on arguments
    if iface_residues_list is not None:
        use_pesto = False
        bound_ref = None
    elif args.apo:
        use_pesto = True
        bound_ref = None
        print("  Apo structure detected - using PeSTo method")
    else:
        use_pesto = False
        bound_ref = args.structure
        print("  Holo structure detected - using 5Å method")

    binding_idx, clusters = derive_binding_indices_direct(
        args.structure, bound_ref, args.bound_chain, args.bound_ligand,
        use_pesto, iface_residues_list, args.plddt_thresh)

    # ── Reference embeddings ──
    p = PDBParser(QUIET=True)
    s = p.get_structure('ref', args.structure)
    chain = list(list(s[0].get_chains()))[0]
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    ref_seq = ''.join(seq1(r.resname) for r in residues)

    print(f"Computing reference embeddings (len={len(ref_seq)})...", flush=True)
    ref_emb, ref_contacts = get_outputs(model, batch_converter, ref_seq, 'reference')
    ref_cs = ref_contact_score(ref_contacts, clusters)
    print(f"  Reference contact score: {ref_cs:.4f}", flush=True)

    # ── Candidates ──
    if not os.path.exists(args.candidates_tsv):
        sys.exit(f"ERROR: candidates TSV not found at {args.candidates_tsv}")
    df = pd.read_csv(args.candidates_tsv, sep='\t')
    rows = []

    for i, row in df.iterrows():
        name = row['candidate']
        seq = str(row['sequence_chain1'])
        print(f"[{i+1:2d}/{len(df)}] {name} (len={len(seq)})...", end=' ', flush=True)

        try:
            emb, contacts = get_outputs(model, batch_converter, seq, name)

            cluster_matches = [match_clusters(ref_emb, emb, contacts, c)
                               for c in clusters]
            all_matched = list(dict.fromkeys(
                [idx for m in cluster_matches for idx in m]))

            cs = contact_score_matched(contacts, cluster_matches)
            es = chem_sim_matched(ref_emb, emb, binding_idx, all_matched)

            es_str = f"{es:.4f}" if not np.isnan(es) else 'N/A'
            print(f"contact={cs:.4f}  chem_sim={es_str}", flush=True)
        except Exception as e:
            cs, es = float('nan'), float('nan')
            print(f"ERROR: {e}", flush=True)

        rows.append({
            'candidate': name,
            'seq_len': len(seq),
            'contact_score': round(cs, 4),
            'chem_sim': round(es, 4) if not np.isnan(es) else None,
        })

    results_df = pd.DataFrame(rows)

    cs_thresh = ref_cs * 0.50
    es_thresh = 0.80

    out_dir = os.path.dirname(out_path)
    if out_dir:  # Only create directory if path contains one
        os.makedirs(out_dir, exist_ok=True)
    results_df = results_df.drop('contact_score', axis=1)
    results_df.to_csv(out_path, sep='\t', index=False)

    print(f"\n{'='*70}")
    print(f"Reference contact score: {ref_cs:.4f}")
    print(f"Contact threshold (50%): {cs_thresh:.4f}")
    print(f"Chemistry sim threshold: {es_thresh:.2f}")
    print(f"{'='*70}")

    # Compute pass/fail on the fly for console output only

    chem_pass = results_df.apply(
        lambda r: True if r['chem_sim'] is None else r['chem_sim'] >= es_thresh, axis=1)
    pass_filter =  chem_pass

    passed_count = pass_filter.sum()
    filtered_count = len(results_df) - passed_count
    print(f"\nPASS: {passed_count}  |  FILTERED: {filtered_count}")

    print("\n── FILTERED OUT ──")
    for i, r in results_df.iterrows():
        if not pass_filter.iloc[i]:
            reasons = []
            if not chem_pass.iloc[i] and r['chem_sim'] is not None:
                reasons.append(f"chem_sim={r['chem_sim']:.4f}<{es_thresh:.2f}")
            print(f"  {r['candidate']:50s}  {', '.join(reasons)}")

    print(f"\nResults saved to {out_path}")

if __name__ == '__main__':
    main()
