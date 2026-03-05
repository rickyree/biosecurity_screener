"""
get_interface_res.py — identify interface residues on an apo receptor using two methods.

Workflow 1  (<5Å):
  Superpose apo onto holo, transform binding partner coords into apo frame,
  flag every apo residue with any atom within 5Å of the partner.

Workflow 2  (PeSTo + P2Rank union):
  Clean apo with pdb-cleaner, run PeSTo and P2Rank via the Amina CLI,
  take the union of PeSTo interface residues and P2Rank pocket-1 residues.

Usage:
  python get_interface_res.py --apo APO.pdb --holo HOLO.pdb [options]

Partner detection (auto, same logic as dmasif_workflow.py):
  1 protein chain in holo  → ligand interface; --bound-ligand required
  2 protein chains in holo → protein-protein interface
  >2 protein chains        → error

Output:
  Prints residue lists for both methods.
  Saves to --out TSV (default: results/interface_residues.tsv).
"""

import os, sys, argparse, subprocess, tempfile, json
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.SeqUtils import seq1
from scipy.spatial import KDTree
import pandas as pd


AMINA = os.path.join(os.path.dirname(sys.executable), 'amina')


# ── PDB helpers ───────────────────────────────────────────────────────────────

def protein_chains(structure):
    chains = []
    for model in structure:
        for chain in model:
            has_protein = any(
                res.id[0] == ' ' and 'CA' in res
                for res in chain
            )
            if has_protein:
                chains.append(chain.id)
        break
    return chains


def ca_residues(structure, chain_id=None):
    """Return {resnum: residue} for standard residues with Cα."""
    out = {}
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0] == ' ' and 'CA' in res:
                    out[res.id[1]] = res
        break
    return out


def resolve_holo_partner(holo_struct, bound_ligand=None):
    """Return (receptor_chain, partner_type, partner_info).

    partner_type: 'ligand' or 'protein'
    partner_info: ligand resname or partner chain id
    """
    chains = protein_chains(holo_struct)

    if len(chains) > 2:
        sys.exit(
            f"ERROR: holo structure has {len(chains)} protein chains. "
            "Only 1 (ligand-bound) or 2 (protein-protein) are supported."
        )

    if len(chains) == 1:
        receptor = chains[0]
        if not bound_ligand:
            sys.exit(
                "ERROR: holo has one protein chain — specify the ligand "
                "residue name with --bound-ligand."
            )
        return receptor, 'ligand', bound_ligand

    # Two protein chains — receptor is chain A if present, else first
    receptor = 'A' if 'A' in chains else chains[0]
    partner  = next(c for c in chains if c != receptor)
    return receptor, 'protein', partner


def get_partner_coords(holo_struct, receptor_chain, partner_type, partner_info):
    """Return (N,3) array of partner heavy-atom coords."""
    coords = []
    for model in holo_struct:
        for chain in model:
            if partner_type == 'ligand':
                if chain.id != receptor_chain:
                    continue
                for res in chain:
                    if res.get_resname().strip() == partner_info:
                        coords.extend(a.coord for a in res)
            else:  # protein
                if chain.id != partner_info:
                    continue
                for res in chain:
                    if res.id[0] == ' ':
                        coords.extend(
                            a.coord for a in res
                            if (a.element or a.name[0]).strip().upper() != 'H'
                        )
        break
    if not coords:
        sys.exit(
            f"ERROR: no atoms found for partner "
            f"({'ligand ' + partner_info if partner_type=='ligand' else 'chain ' + partner_info})."
        )
    return np.array(coords)


# ── Workflow 1: <5Å ──────────────────────────────────────────────────────────

def workflow_5a(apo_path, holo_path, bound_ligand=None, cutoff=5.0):
    """Return sorted list of apo residue numbers within cutoff of holo partner."""
    parser = PDBParser(QUIET=True)
    holo   = parser.get_structure('holo', holo_path)
    apo    = parser.get_structure('apo',  apo_path)

    receptor_chain, partner_type, partner_info = resolve_holo_partner(holo, bound_ligand)
    partner_coords = get_partner_coords(holo, receptor_chain, partner_type, partner_info)
    print(f"  Partner: {partner_type} '{partner_info}'  ({len(partner_coords)} heavy atoms)")

    # Superpose apo onto holo receptor by shared residue number
    holo_res = ca_residues(holo, receptor_chain)
    apo_res  = ca_residues(apo)
    common   = sorted(set(holo_res) & set(apo_res))
    if len(common) < 10:
        sys.exit(f"ERROR: only {len(common)} shared Cα residues — cannot superpose.")

    holo_ca = [holo_res[i]['CA'] for i in common]
    apo_ca  = [apo_res[i]['CA']  for i in common]

    sup = Superimposer()
    sup.set_atoms(holo_ca, apo_ca)   # fixed=holo, moving=apo
    sup.apply(list(apo.get_atoms()))
    print(f"  Superposition: {len(common)} Cα pairs, RMSD = {sup.rms:.3f} Å")

    # Partner is already in holo frame; apo is now also in holo frame
    tree = KDTree(partner_coords)
    interface = []
    for resnum, res in sorted(apo_res.items()):
        for atom in res:
            if tree.query_ball_point(atom.coord, r=cutoff):
                interface.append(resnum)
                break

    return sorted(interface)


# ── Workflow 2: PeSTo + P2Rank union ─────────────────────────────────────────

def _amina_run(args, background=False):
    """Run an amina CLI command. Returns job-id if background, else None."""
    cmd = [AMINA, 'run'] + args
    if background:
        cmd.append('--background')
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr
    if result.returncode != 0 and 'Job submitted' not in out:
        print(f"  [amina] WARNING: {out.strip()[:300]}", file=sys.stderr)
        return None
    # Extract job-id from output
    for line in out.splitlines():
        if 'Job submitted:' in line:
            return line.split('Job submitted:')[1].strip()
    return None


def _amina_wait_download(job_id, out_dir):
    subprocess.run([AMINA, 'jobs', 'wait', job_id], check=True)
    subprocess.run([AMINA, 'jobs', 'download', job_id, '-o', out_dir], check=True)


def workflow_union(apo_path, holo_path, bound_ligand=None,
                   pesto_threshold=0.5, work_dir='results/interface_tmp'):
    """Return sorted list of residue numbers from PeSTo ∪ P2Rank pocket-1."""
    os.makedirs(work_dir, exist_ok=True)
    parser = PDBParser(QUIET=True)
    holo   = parser.get_structure('holo', holo_path)
    _, partner_type, _ = resolve_holo_partner(holo, bound_ligand)

    # Step 1: clean apo
    print("  [union] Cleaning apo with pdb-cleaner …", flush=True)
    clean_dir  = os.path.join(work_dir, 'cleaned')
    apo_stem   = os.path.splitext(os.path.basename(apo_path))[0]
    clean_path = os.path.join(clean_dir, f'{apo_stem}_cleaned.pdb')

    if not os.path.exists(clean_path):
        result = subprocess.run(
            [AMINA, 'run', 'pdb-cleaner', '--pdb', apo_path,
             '-o', clean_dir, '-j', apo_stem],
            capture_output=True, text=True
        )
        if not os.path.exists(clean_path):
            print(f"  [union] pdb-cleaner failed; using original apo.", file=sys.stderr)
            clean_path = apo_path
    else:
        print(f"  [union] Using cached cleaned apo: {clean_path}")

    # Step 2: submit PeSTo and P2Rank in parallel
    pesto_iface = 'Protein-Ligand-Interface' if partner_type == 'ligand' \
                  else 'Protein-Protein-Interface'
    print(f"  [union] Submitting PeSTo ({pesto_iface}) …", flush=True)
    pesto_dir  = os.path.join(work_dir, 'pesto')
    p2rank_dir = os.path.join(work_dir, 'p2rank')
    os.makedirs(pesto_dir,  exist_ok=True)
    os.makedirs(p2rank_dir, exist_ok=True)

    pesto_job  = _amina_run(['pesto', '--pdb', clean_path,
                              '--threshold', str(pesto_threshold),
                              '-o', pesto_dir, '-j', apo_stem], background=True)
    print(f"  [union] Submitting P2Rank …", flush=True)
    p2rank_job = _amina_run(['p2rank', '--pdb', clean_path,
                              '-o', p2rank_dir, '-j', apo_stem], background=True)

    # Step 3: wait and download
    for job_id, label, out_dir in [
        (pesto_job,  'PeSTo',  pesto_dir),
        (p2rank_job, 'P2Rank', p2rank_dir),
    ]:
        if job_id:
            print(f"  [union] Waiting for {label} (job {job_id[:8]}) …", flush=True)
            _amina_wait_download(job_id, out_dir)

    # Step 4: parse PeSTo — the relevant interface type CSV
    iface_key   = 'Protein-Ligand-Interface' if partner_type == 'ligand' \
                  else 'Protein-Protein-Interface'
    pesto_csv   = os.path.join(pesto_dir,
                               f'pesto_{apo_stem}_bfactor_{iface_key}_residues.csv')
    if not os.path.exists(pesto_csv):
        # fall back to any residues CSV
        csvs = [f for f in os.listdir(pesto_dir) if f.endswith('_residues.csv')]
        pesto_csv = os.path.join(pesto_dir, csvs[0]) if csvs else None

    pesto_res = set()
    if pesto_csv and os.path.exists(pesto_csv):
        df = pd.read_csv(pesto_csv)
        pesto_res = set(df['residue_number'].tolist())
        print(f"  [union] PeSTo ({iface_key}): {len(pesto_res)} residues")
    else:
        print(f"  [union] WARNING: PeSTo CSV not found.", file=sys.stderr)

    # Step 5: parse P2Rank — pocket 1 residues only
    p2rank_csv = os.path.join(p2rank_dir, f'p2rank_{apo_stem}_residues.csv')
    p2rank_res = set()
    if os.path.exists(p2rank_csv):
        df = pd.read_csv(p2rank_csv, skipinitialspace=True)
        p2rank_res = set(df[df['pocket'] == 1]['residue_label'].tolist())
        print(f"  [union] P2Rank pocket-1: {len(p2rank_res)} residues")
    else:
        print(f"  [union] WARNING: P2Rank residues CSV not found.", file=sys.stderr)

    union = pesto_res | p2rank_res
    return sorted(union)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--apo',          required=True, metavar='PDB',
                    help='Apo receptor PDB (protein A).')
    ap.add_argument('--holo',         required=True, metavar='PDB',
                    help='Holo/complex PDB defining the binding site.')
    ap.add_argument('--bound-ligand', default=None, metavar='RES',
                    help='Ligand residue name in holo (required for single-chain holo).')
    ap.add_argument('--cutoff',       type=float, default=5.0,
                    help='Distance cutoff in Å for workflow 1 (default: 5.0).')
    ap.add_argument('--pesto-threshold', type=float, default=0.5,
                    help='PeSTo binding probability threshold (default: 0.5).')
    ap.add_argument('--work-dir',     default='results/interface_tmp',
                    help='Working directory for Amina outputs.')
    ap.add_argument('--out',          default='results/interface_residues.tsv',
                    help='Output TSV path.')
    ap.add_argument('--skip-5a',      action='store_true',
                    help='Skip workflow 1 (<5Å).')
    ap.add_argument('--skip-union',   action='store_true',
                    help='Skip workflow 2 (PeSTo+P2Rank union).')
    args = ap.parse_args()

    apo_stem = os.path.splitext(os.path.basename(args.apo))[0]

    res_5a    = []
    res_union = []

    if not args.skip_5a:
        print("\n── Workflow 1: <5Å from holo binding partner ──────────────────")
        res_5a = workflow_5a(args.apo, args.holo, args.bound_ligand, args.cutoff)
        print(f"  Result ({len(res_5a)} residues): {res_5a}")

    if not args.skip_union:
        print("\n── Workflow 2: PeSTo + P2Rank union ───────────────────────────")
        res_union = workflow_union(args.apo, args.holo, args.bound_ligand,
                                   args.pesto_threshold, args.work_dir)
        print(f"  Result ({len(res_union)} residues): {res_union}")

    # Summary
    if res_5a and res_union:
        overlap = sorted(set(res_5a) & set(res_union))
        print(f"\n── Overlap between methods: {len(overlap)} residues ──────────────")
        print(f"  {overlap}")

    # Save TSV
    os.makedirs(os.path.dirname(args.out) if os.path.dirname(args.out) else '.', exist_ok=True)
    rows = []
    all_res = sorted(set(res_5a) | set(res_union))
    for r in all_res:
        rows.append({
            'residue_number': r,
            'in_5a':          r in res_5a,
            'in_union':       r in res_union,
            'in_both':        r in res_5a and r in res_union,
        })
    pd.DataFrame(rows).to_csv(args.out, sep='\t', index=False)
    print(f"\nSaved → {args.out}")


if __name__ == '__main__':
    main()
