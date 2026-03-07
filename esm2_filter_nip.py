#!/usr/bin/env python
"""
ESM2 filtering workflow for nipah sequences using explicit interface residues.
Based on esm2_prefilter.py but uses explicit residue indices instead of dMaSIF cache.
"""

import sys, os, hashlib
import numpy as np
import pandas as pd
import torch
import esm as esmlib
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ESM2_CACHE_DIR = 'results/esm2_embed_cache'

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

def match_clusters(ref_emb, cand_emb, cand_contacts, ref_cluster_idx, contact_thresh=0.1):
    """For each reference cluster index, find the nearest-neighbour position
    in the candidate via cosine similarity in embedding space."""
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
    their best match within the already-matched, compacted candidate patch."""
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

def main():
    # Interface residues from <5A analysis (1-indexed from PDB)
    interface_residues = [51, 52, 53, 116, 199, 200, 212, 213, 269, 299, 300, 301, 302, 303, 312,
                         315, 316, 317, 318, 341, 342, 343, 344, 366, 367, 368, 369, 370, 390, 391, 392, 394, 399]

    print("Loading ESM2 650M...", flush=True)
    model, alphabet = esmlib.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model  = model.to(device)
    model.eval()
    print(f"  device: {device}", flush=True)

    # ── Reference setup ──
    ref_pdb = 'binding/2vsm.pdb'
    print(f"Reference: {ref_pdb} chain A", flush=True)

    parser = PDBParser(QUIET=True)
    s = parser.get_structure('ref', ref_pdb)
    chain = list(s[0].get_chains())[0]  # Chain A
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    ref_seq = ''.join(seq1(r.resname) for r in residues)

    print(f"Reference sequence length: {len(ref_seq)}")

    # Convert PDB residue numbers (1-indexed) to sequence indices (0-indexed)
    binding_idx = [r-1 for r in interface_residues if r <= len(ref_seq)]
    print(f"Interface residues (0-indexed): {len(binding_idx)} residues")
    print(f"Interface indices: {binding_idx}")

    # Create clusters - for simplicity, treat as one large cluster
    clusters = [binding_idx]

    print(f"Computing reference embeddings (len={len(ref_seq)})...", flush=True)
    ref_emb, ref_contacts = get_outputs(model, batch_converter, ref_seq, 'reference')
    ref_cs = ref_contact_score(ref_contacts, clusters)
    print(f"  Reference contact score: {ref_cs:.4f}", flush=True)

    # ── Load and filter candidates ──
    candidates_file = 'results/dmasif/candidate_sequences.tsv'
    print(f"Loading candidates from {candidates_file}")

    df = pd.read_csv(candidates_file, sep='\t')
    # Filter for nip_mut* sequences only
    nip_df = df[df['candidate'].str.startswith('nip_mut')].copy()
    print(f"Found {len(nip_df)} nip_mut* sequences")

    rows = []

    for i, row in nip_df.iterrows():
        name = row['candidate']
        seq  = str(row['sequence_chain1'])
        print(f"[{i+1-nip_df.index[0]+1:2d}/{len(nip_df)}] {name} (len={len(seq)})...", end=' ', flush=True)

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

    output_file = 'results/dmasif/esm2_prefilter_nip.tsv'
    results_df.to_csv(output_file, sep='\t', index=False)

    print(f"\n{'='*70}")
    print(f"Reference contact score: {ref_cs:.4f}")
    print(f"Contact threshold (50%): {cs_thresh:.4f}")
    print(f"Chemistry sim threshold: {es_thresh:.2f}")
    print(f"Interface residues used: {len(binding_idx)}")
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

    print(f"\nResults saved to {output_file}")

if __name__ == '__main__':
    main()
