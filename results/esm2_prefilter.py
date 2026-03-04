"""
ESM2 pre-filter for candidate sequences.
Two filters:
  1. Inter-cluster contact score: ESM2 contact probability between sequence-distant
     binding site clusters (C1:idx 19-68, C2:idx 147-152, C3:idx 432-444)
  2. Chemistry similarity: cosine similarity of ESM2 embeddings at 29 binding site
     positions vs. the reference (1f0l_A chain B)
"""
import sys, os
import numpy as np
import pandas as pd
import torch
import esm as esmlib

# ── Binding site indices (0-based, from 1f0l_A chain B residue→seqidx mapping) ──
BINDING_IDX = [
    19,20,21,22,23,26,29,30,33,34,35,37,41,42,43,44,45,
    51,52,53,54,55,61,64,68,   # cluster 1 (residues 20-69)
    147,152,                    # cluster 2 (residues 148,153)
    432,444                     # cluster 3 (residues 446,458)
]
C1 = [19,20,21,22,23,26,29,30,33,34,35,37,41,42,43,44,45,51,52,53,54,55,61,64,68]
C2 = [147,152]
C3 = [432,444]

MIN_LEN_FOR_C3 = 445   # need at least index 444


def get_outputs(model, batch_converter, seq, label):
    clean = seq.replace('X', 'G').replace('-', 'G')[:1000]  # ESM2 safe length
    data = [(label, clean)]
    _, _, tokens = batch_converter(data)
    device = next(model.parameters()).device
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[33], return_contacts=True)
    emb = out['representations'][33][0, 1:1+len(clean)].cpu().numpy()
    contacts = out['contacts'][0].cpu().numpy()
    return emb, contacts


def contact_score(contacts, c1, c2, c3, L):
    valid_c1 = [i for i in c1 if i < L]
    valid_c2 = [i for i in c2 if i < L]
    valid_c3 = [i for i in c3 if i < L]
    pairs = []
    if valid_c1 and valid_c3:
        pairs.append(contacts[np.ix_(valid_c1, valid_c3)].mean())
    if valid_c2 and valid_c3:
        pairs.append(contacts[np.ix_(valid_c2, valid_c3)].mean())
    if valid_c1 and valid_c2:
        pairs.append(contacts[np.ix_(valid_c1, valid_c2)].mean())
    return float(np.mean(pairs)) if pairs else float('nan')


def chem_sim(emb, ref_emb, binding_idx):
    valid = [i for i in binding_idx if i < len(emb)]
    if len(valid) < 10:
        return float('nan')
    ce = emb[valid]
    re = ref_emb[valid]
    cn = np.linalg.norm(ce, axis=1, keepdims=True)
    rn = np.linalg.norm(re, axis=1, keepdims=True)
    cn[cn == 0] = 1e-8
    rn[rn == 0] = 1e-8
    return float(((ce / cn) * (re / rn)).sum(axis=1).mean())


def main():
    # ── Load model ──
    print("Loading ESM2 650M...", flush=True)
    model, alphabet = esmlib.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    device = 'mps' if torch.backends.mps.is_available() else 'cpu'
    model = model.to(device)
    model.eval()
    print(f"  device: {device}", flush=True)

    # ── Reference sequence from 1f0l_A chain B ──
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1
    p = PDBParser(QUIET=True)
    s = p.get_structure('ref', 'binding/1f0l_A.pdb')
    chain = list(list(s[0].get_chains()))[0]
    residues = list(chain.get_residues())
    ref_seq = ''.join(seq1(r.resname) for r in residues)

    print(f"Computing reference embeddings (len={len(ref_seq)})...", flush=True)
    ref_emb, ref_contacts = get_outputs(model, batch_converter, ref_seq, 'reference')
    ref_cs = contact_score(ref_contacts, C1, C2, C3, len(ref_seq))
    print(f"  Reference contact score: {ref_cs:.4f}", flush=True)

    # ── Candidate sequences ──
    df = pd.read_csv('results/dmasif/candidate_sequences.tsv', sep='\t')
    rows = []

    for i, row in df.iterrows():
        name = row['candidate']
        seq = str(row['sequence_chain1'])
        L = len(seq.replace('X','G').replace('-','G')[:1000])
        is_scaffold = name.startswith('1f0l_A_mut') or name.startswith('1f0l_mut')

        print(f"[{i+1:2d}/{len(df)}] {name} (len={len(seq)})...", end=' ', flush=True)

        try:
            emb, contacts = get_outputs(model, batch_converter, seq, name)
            cs = contact_score(contacts, C1, C2, C3, L)
            es = chem_sim(emb, ref_emb, BINDING_IDX) if is_scaffold else float('nan')
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
            'scaffold': is_scaffold,
        })

    results_df = pd.DataFrame(rows)

    # ── Thresholds ──
    # Contact: flag if < 30% of reference contact score
    cs_thresh = ref_cs * 0.30
    # Chem sim: flag if < 0.80 (80% cosine similarity)
    es_thresh = 0.80

    results_df['contact_pass'] = results_df['contact_score'].apply(
        lambda x: True if np.isnan(x) else x >= cs_thresh)
    results_df['chem_pass'] = results_df.apply(
        lambda r: True if r['chem_sim'] is None else r['chem_sim'] >= es_thresh, axis=1)
    results_df['pass_filter'] = results_df['contact_pass'] & results_df['chem_pass']

    out_path = 'results/dmasif/esm2_prefilter.tsv'
    results_df.to_csv(out_path, sep='\t', index=False)

    print(f"\n{'='*70}")
    print(f"Reference contact score: {ref_cs:.4f}")
    print(f"Contact threshold (30%): {cs_thresh:.4f}")
    print(f"Chemistry sim threshold: {es_thresh:.2f}")
    print(f"{'='*70}")

    # Summary
    filtered = results_df[~results_df['pass_filter']]
    passed = results_df[results_df['pass_filter']]
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
