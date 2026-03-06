"""
ESM2 pre-filter for candidate sequences.
Two filters (both position-free: reference binding site embeddings are matched
to candidate positions via nearest-neighbour search in embedding space):
  1. Inter-cluster contact score: ESM2 contact probability between candidate
     positions that best match the reference binding site cluster embeddings.
  2. Chemistry similarity: mean cosine similarity of reference binding site
     embeddings to their nearest neighbours in the candidate patch.

Binding site indices and clusters are derived from the dMaSIF interface cache
(A_interface_cache_{apo_stem}_{bound_stem}.npz), which must be built first by
running dmasif_workflow.py with the same --apo-ref and --bound-ref arguments.
"""
import sys, os, argparse, hashlib
import numpy as np
import pandas as pd
import torch
import esm as esmlib
from scipy.spatial import KDTree

ESM2_CACHE_DIR = 'results/esm2_embed_cache'


# ── ESM2 inference (with sequence-keyed cache) ────────────────────────────────

def seq_hash(seq):
    """SHA256 hash of the cleaned sequence — used as cache key."""
    clean = seq.replace('X', 'G').replace('-', 'G')[:1000]
    return hashlib.sha256(clean.encode()).hexdigest()


def get_outputs(model, batch_converter, seq, label):
    clean     = seq.replace('X', 'G').replace('-', 'G')[:1000]
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
    emb      = out['representations'][33][0, 1:1+len(clean)].cpu().numpy()
    contacts = out['contacts'][0].cpu().numpy()

    os.makedirs(ESM2_CACHE_DIR, exist_ok=True)
    np.savez(cache_path, emb=emb, contacts=contacts)
    return emb, contacts


# ── SGK binding site derivation ───────────────────────────────────────────────

def derive_binding_indices(cache_path, apo_pdb, gap=10):
    """Load a dMaSIF interface cache, map iface_coords to apo Cα sequence
    indices, and split into clusters by sequence gaps >= `gap`.

    Returns (binding_idx, clusters) where clusters is a list of index lists.
    """
    from Bio.PDB import PDBParser

    if not os.path.exists(cache_path):
        sys.exit(f"ERROR: dMaSIF cache not found at {cache_path}. "
                 f"Run dmasif_workflow.py with matching --apo-ref and --bound-ref first.")
    if not os.path.exists(apo_pdb):
        sys.exit(f"ERROR: apo PDB not found at {apo_pdb}.")

    cache        = np.load(cache_path)
    iface_coords = cache['A_iface_coords']   # (N, 3) coords in apo frame

    p        = PDBParser(QUIET=True)
    apo_s    = p.get_structure('APO', apo_pdb)
    residues = [res for model in apo_s for chain in model
                for res in chain if res.id[0] == ' ' and 'CA' in res]
    ca_coords = np.array([res['CA'].coord for res in residues])

    tree    = KDTree(ca_coords)
    _, idx  = tree.query(iface_coords)
    seq_idx = sorted(set(idx.tolist()))

    clusters, current = [], [seq_idx[0]]
    for i in range(1, len(seq_idx)):
        if seq_idx[i] - seq_idx[i-1] >= gap:
            clusters.append(current)
            current = []
        current.append(seq_idx[i])
    clusters.append(current)

    print(f"  Binding site: {len(seq_idx)} residues, "
          f"{len(clusters)} cluster(s) starting at: {[c[0] for c in clusters]}")

    return seq_idx, clusters


# ── Matching and scoring ──────────────────────────────────────────────────────

def match_clusters(ref_emb, cand_emb, cand_contacts, ref_cluster_idx,
                   contact_thresh=0.1):
    """For each reference cluster index, find the nearest-neighbour position
    in the candidate via cosine similarity in embedding space.

    Compaction uses the candidate's ESM2 contact map rather than sequence
    distance: a matched position is kept only if its mean contact probability
    to the other matched positions meets `contact_thresh`. This correctly
    handles binding sites where relevant residues are far apart in sequence
    but spatially proximal in structure (e.g. P-loop, hinge, DFG in kinases).
    """
    cand_norm = cand_emb / (np.linalg.norm(cand_emb, axis=1, keepdims=True) + 1e-8)

    matched = []
    for ri in ref_cluster_idx:
        if ri >= len(ref_emb):
            continue
        rv      = ref_emb[ri]
        rv_norm = rv / (np.linalg.norm(rv) + 1e-8)
        matched.append(int(np.argmax(cand_norm @ rv_norm)))

    if not matched:
        return []

    matched = list(dict.fromkeys(matched))   # deduplicate

    if len(matched) == 1:
        return matched

    # Contact-based compaction: keep positions mutually proximal in structure.
    # A position survives if its mean ESM2 contact probability to the other
    # matched positions is >= contact_thresh.
    L = cand_contacts.shape[0]
    valid = [i for i in matched if i < L]
    if len(valid) < 2:
        return valid

    sub = cand_contacts[np.ix_(valid, valid)]   # (M, M) contact submatrix
    np.fill_diagonal(sub, 0.0)
    mean_contact = sub.mean(axis=1)             # mean contact to peers
    compact = [valid[i] for i in range(len(valid))
               if mean_contact[i] >= contact_thresh]

    return compact if compact else valid         # fallback: keep all if none pass


def contact_score_matched(contacts, cluster_matches):
    """Inter-cluster contact score across all pairs of matched cluster sets."""
    L     = contacts.shape[0]
    valid = [[i for i in m if i < L] for m in cluster_matches]
    pairs = []
    for a in range(len(valid)):
        for b in range(a + 1, len(valid)):
            if valid[a] and valid[b]:
                pairs.append(contacts[np.ix_(valid[a], valid[b])].mean())
    return float(np.mean(pairs)) if pairs else float('nan')


def chem_sim_matched(ref_emb, cand_emb, ref_binding_idx, matched_cand_idx):
    """Mean cosine similarity between reference binding site embeddings and
    their best match within the already-matched, compacted candidate patch.

    Restricting search to `matched_cand_idx` removes the large-protein bias
    from searching all positions.
    """
    if len(matched_cand_idx) < 3:
        return float('nan')

    patch      = cand_emb[matched_cand_idx]
    patch_norm = patch / (np.linalg.norm(patch, axis=1, keepdims=True) + 1e-8)

    sims = []
    for ri in ref_binding_idx:
        if ri >= len(ref_emb):
            continue
        rv      = ref_emb[ri]
        rv_norm = rv / (np.linalg.norm(rv) + 1e-8)
        sims.append(float((patch_norm @ rv_norm).max()))

    if len(sims) < 3:
        return float('nan')
    return float(np.mean(sims))


def ref_contact_score(ref_contacts, clusters):
    """Contact score on the reference itself using its own cluster indices."""
    L     = ref_contacts.shape[0]
    valid = [[i for i in c if i < L] for c in clusters]
    pairs = []
    for a in range(len(valid)):
        for b in range(a + 1, len(valid)):
            if valid[a] and valid[b]:
                pairs.append(ref_contacts[np.ix_(valid[a], valid[b])].mean())
    return float(np.mean(pairs)) if pairs else float('nan')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--apo-ref', required=True, metavar='PDB',
                    help='Apo receptor PDB — must match the --apo-ref used in dmasif_workflow.py.')
    ap.add_argument('--bound-ref', default=None, metavar='PDB',
                    help='Bound reference PDB — must match the --bound-ref used in dmasif_workflow.py. '
                         'Not required when --union or --iface-residues was used.')
    ap.add_argument('--candidates-tsv', required=True, metavar='TSV',
                    help='TSV file of candidate sequences with columns '
                         '"candidate" and "sequence_chain1".')
    ap.add_argument('--out', default=None, metavar='TSV',
                    help='Output TSV path (default: results/dmasif/esm2_prefilter_{apo_stem}[_union|_explicit].tsv)')
    ap.add_argument('--union', action='store_true',
                    help='Use the union interface cache (built with --union in dmasif_workflow.py).')
    ap.add_argument('--iface-residues', action='store_true',
                    help='Use the explicit-residues cache (built with --iface-residues in dmasif_workflow.py).')
    args = ap.parse_args()

    apo_stem   = os.path.splitext(os.path.basename(args.apo_ref))[0]
    bound_stem = os.path.splitext(os.path.basename(args.bound_ref))[0] if args.bound_ref else 'noref'
    if args.iface_residues:
        suffix = '_explicit'
    elif args.union:
        suffix = '_union'
    else:
        suffix = ''
    cache_path = f'results/dmasif/A_interface_cache_{apo_stem}_{bound_stem}{suffix}.npz'
    out_path   = args.out or f'results/dmasif/esm2_prefilter_{apo_stem}{suffix}.tsv'

    print("Loading ESM2 650M...", flush=True)
    model, alphabet = esmlib.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model  = model.to(device)
    model.eval()
    print(f"  device: {device}", flush=True)

    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1
    p = PDBParser(QUIET=True)

    # ── Reference setup ──
    print(f"Reference: {args.apo_ref}", flush=True)
    binding_idx, clusters = derive_binding_indices(cache_path, args.apo_ref)
    ref_pdb = args.apo_ref

    s        = p.get_structure('ref', ref_pdb)
    chain    = list(list(s[0].get_chains()))[0]
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    ref_seq  = ''.join(seq1(r.resname) for r in residues)

    print(f"Computing reference embeddings (len={len(ref_seq)})...", flush=True)
    ref_emb, ref_contacts = get_outputs(model, batch_converter, ref_seq, 'reference')
    ref_cs = ref_contact_score(ref_contacts, clusters)
    print(f"  Reference contact score: {ref_cs:.4f}", flush=True)

    # ── Candidates ──
    if not os.path.exists(args.candidates_tsv):
        sys.exit(f"ERROR: candidates TSV not found at {args.candidates_tsv}")
    df   = pd.read_csv(args.candidates_tsv, sep='\t')
    rows = []

    for i, row in df.iterrows():
        name = row['candidate']
        seq  = str(row['sequence_chain1'])
        print(f"[{i+1:2d}/{len(df)}] {name} (len={len(seq)})...", end=' ', flush=True)

        try:
            emb, contacts = get_outputs(model, batch_converter, seq, name)

            cluster_matches = [match_clusters(ref_emb, emb, contacts, c)
                               for c in clusters]
            all_matched     = list(dict.fromkeys(
                [idx for m in cluster_matches for idx in m]))

            cs = contact_score_matched(contacts, cluster_matches)
            es = chem_sim_matched(ref_emb, emb, binding_idx, all_matched)

            es_str = f"{es:.4f}" if not np.isnan(es) else 'N/A'
            print(f"contact={cs:.4f}  chem_sim={es_str}", flush=True)
        except Exception as e:
            cs, es = float('nan'), float('nan')
            print(f"ERROR: {e}", flush=True)

        rows.append({
            'candidate'    : name,
            'seq_len'      : len(seq),
            'contact_score': round(cs, 4),
            'chem_sim'     : round(es, 4) if not np.isnan(es) else None,
        })

    results_df = pd.DataFrame(rows)

    cs_thresh = ref_cs * 0.50
    es_thresh = 0.80

    results_df['contact_pass'] = results_df['contact_score'].apply(
        lambda x: True if np.isnan(x) else x >= cs_thresh)
    results_df['chem_pass'] = results_df.apply(
        lambda r: True if r['chem_sim'] is None else r['chem_sim'] >= es_thresh, axis=1)
    results_df['pass_filter'] = results_df['contact_pass'] & results_df['chem_pass']

    results_df.to_csv(out_path, sep='\t', index=False)

    print(f"\n{'='*70}")
    print(f"Reference contact score: {ref_cs:.4f}")
    print(f"Contact threshold (50%): {cs_thresh:.4f}")
    print(f"Chemistry sim threshold: {es_thresh:.2f}")
    print(f"{'='*70}")

    filtered = results_df[~results_df['pass_filter']]
    passed   = results_df[results_df['pass_filter']]
    print(f"\nPASS: {len(passed)}  |  FILTERED: {len(filtered)}")

    print("\n── FILTERED OUT ──")
    for _, r in filtered.iterrows():
        reasons = []
        if not r['contact_pass']:
            reasons.append(f"contact={r['contact_score']:.4f}<{cs_thresh:.4f}")
        if not r['chem_pass'] and r['chem_sim'] is not None:
            reasons.append(f"chem_sim={r['chem_sim']:.4f}<{es_thresh:.2f}")
        print(f"  {r['candidate']:50s}  {', '.join(reasons)}")

    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
