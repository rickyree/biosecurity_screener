#!/usr/bin/env python3
"""
ProteinMPNN ligand-site design via Amina CLI.

For a given PDB file:
  1. Identifies protein residues within --cutoff Å of a specified ligand.
  2. Submits N_OUTSIDE ProteinMPNN jobs with those residues fixed (binding site
     conserved, scaffold varied) across a temperature spectrum.
  3. Submits N_CONSERVED ProteinMPNN jobs with no fixed residues (everything
     free to mutate) across a temperature spectrum.
  4. Waits for all jobs, downloads results, parses FASTA outputs.
  5. Appends designed sequences to a TSV file.

Usage:
    python proteinmpnn_ligand_design.py \\
        --pdb binding/1f0l_A.pdb \\
        --chain B \\
        --ligand APU \\
        --prefix 1f0l_A \\
        --n-outside 15 \\
        --n-conserved 10 \\
        --cutoff 5.0 \\
        --output-tsv results/dmasif/candidate_sequences.tsv
"""

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


# ── Interface detection ───────────────────────────────────────────────────────

def find_ligand_interface(pdb_path, chain_id, ligand_name, cutoff=5.0):
    """Return sorted list of residue numbers within cutoff Å of the ligand."""
    from scipy.spatial import KDTree
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('s', pdb_path)
    chain  = struct[0][chain_id]

    ligand_coords = np.array([
        a.coord for r in chain for a in r
        if r.get_resname().strip() == ligand_name
    ])
    if len(ligand_coords) == 0:
        sys.exit(f"ERROR: ligand '{ligand_name}' not found in chain {chain_id} of {pdb_path}")

    tree = KDTree(ligand_coords)
    iface = []
    for res in chain:
        if res.id[0] != ' ':
            continue
        for atom in res:
            if tree.query_ball_point(atom.coord, r=cutoff):
                iface.append(res.id[1])
                break
    return sorted(iface)


def pdb_to_sequential(pdb_path, chain_id, pdb_resnums):
    """Map PDB residue numbers to 1-based sequential positions (as ProteinMPNN expects)."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('s', pdb_path)
    ordered = [r.id[1] for r in struct[0][chain_id]
               if r.id[0] == ' ' and 'CA' in r]
    resnum_to_seq = {rn: i + 1 for i, rn in enumerate(ordered)}
    mapped = [resnum_to_seq[r] for r in pdb_resnums if r in resnum_to_seq]
    missing = [r for r in pdb_resnums if r not in resnum_to_seq]
    if missing:
        print(f"  WARNING: residues not found in chain {chain_id}: {missing}")
    return mapped


def find_protein_interface(pdb_path, chain_id, partner_chain_id, cutoff=5.0):
    """Return sorted list of receptor residue numbers within cutoff Å of partner chain."""
    from scipy.spatial import KDTree
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('s', pdb_path)
    model  = struct[0]

    if partner_chain_id not in [c.id for c in model]:
        sys.exit(f"ERROR: partner chain '{partner_chain_id}' not found in {pdb_path}")

    partner_coords = np.array([
        a.coord for res in model[partner_chain_id]
        if res.id[0] == ' '
        for a in res
        if (a.element or a.name[0]).strip().upper() != 'H'
    ])
    if len(partner_coords) == 0:
        sys.exit(f"ERROR: no heavy atoms in partner chain '{partner_chain_id}' of {pdb_path}")

    tree  = KDTree(partner_coords)
    iface = []
    for res in model[chain_id]:
        if res.id[0] != ' ':
            continue
        for atom in res:
            if tree.query_ball_point(atom.coord, r=cutoff):
                iface.append(res.id[1])
                break
    return sorted(iface)


# ── Job submission ─────────────────────────────────────────────────────────────

def submit_job(pdb_path, chain_id, fixed_str, temperature, job_name):
    """Submit a background ProteinMPNN job and return the job ID."""
    cmd = [
        'amina', 'run', 'proteinmpnn',
        '--pdb', str(pdb_path),
        '--chains', chain_id,
        '-n', '1',
        '--temperature', str(round(temperature, 2)),
        '--job-name', job_name,
        '--background',
    ]
    if fixed_str:
        cmd += ['--fixed', fixed_str]

    result = subprocess.run(cmd, capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if 'Job submitted:' in line:
            return line.split('Job submitted:')[1].strip()

    print(f"  WARNING: could not parse job ID for {job_name}")
    print(result.stdout[-400:])
    return None


def wait_for_jobs(job_ids, poll_interval=10):
    """Block until all job IDs complete."""
    result = subprocess.run(
        ['amina', 'jobs', 'wait'] + job_ids + ['--poll-interval', str(poll_interval)],
        capture_output=True, text=True
    )
    print(result.stdout[-1000:])
    if result.returncode != 0:
        print(result.stderr[-500:])


def download_job(job_id, out_dir):
    """Download job artifacts to out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ['amina', 'jobs', 'download', job_id, '-o', str(out_dir)],
        capture_output=True
    )


# ── Sequence parsing ──────────────────────────────────────────────────────────

def parse_fasta(job_dir):
    """Return the designed sequence from a ProteinMPNN job directory."""
    fastas = list(job_dir.rglob('*.fasta')) + list(job_dir.rglob('*.fa'))
    if not fastas:
        return None
    recs = list(SeqIO.parse(str(fastas[0]), 'fasta'))
    designed = [r for r in recs
                if 'original' not in r.id.lower() and 'native' not in r.id.lower()]
    seq = str(designed[0].seq) if designed else str(recs[0].seq)
    return re.sub(r'[^A-Z]', '', seq)


# ── TSV update ────────────────────────────────────────────────────────────────

def append_to_tsv(tsv_path, sequences):
    """Append {name: sequence} dict to the candidate TSV, creating if needed."""
    tsv_path = Path(tsv_path)
    if tsv_path.exists():
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            header = reader.fieldnames or ['candidate', 'sequence_chain1']
            rows = list(reader)
    else:
        header = ['candidate', 'sequence_chain1']
        rows = []

    for name, seq in sequences.items():
        row = {h: '' for h in header}
        row['candidate'] = name
        row['sequence_chain1'] = seq
        rows.append(row)

    with open(tsv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pdb',           required=True,  help='Input PDB file')
    ap.add_argument('--chain',         required=True,  help='Chain ID of the receptor to redesign')
    ap.add_argument('--prefix',        required=True,  help='Candidate name prefix (e.g. 1f0l_A)')

    iface_grp = ap.add_mutually_exclusive_group(required=True)
    iface_grp.add_argument('--ligand',         metavar='RES',
                            help='Ligand residue name to define binding site (e.g. APU, NAD)')
    iface_grp.add_argument('--partner-chain',  metavar='ID',
                            help='Chain ID of the protein partner to define the interface (e.g. B)')
    iface_grp.add_argument('--iface-residues', metavar='NUMS',
                            help='Comma-separated residue numbers to use as the binding site '
                                 '(e.g. 32,64,65,68). Bypasses proximity detection; '
                                 'useful when --pdb is an apo and residues were identified from a holo.')
    ap.add_argument('--n-outside',    type=int, default=15,
                    help='Number of mutoutside sequences (default: 15)')
    ap.add_argument('--n-conserved',  type=int, default=10,
                    help='Number of mutconserved sequences (default: 10)')
    ap.add_argument('--cutoff',       type=float, default=5.0,
                    help='Distance cutoff in Å for binding site (default: 5.0)')
    ap.add_argument('--temp-min',     type=float, default=0.1,
                    help='Minimum sampling temperature (default: 0.1)')
    ap.add_argument('--temp-max',     type=float, default=1.5,
                    help='Maximum sampling temperature (default: 1.5)')
    ap.add_argument('--output-tsv',   default='results/dmasif/candidate_sequences.tsv',
                    help='TSV file to append sequences to')
    ap.add_argument('--output-dir',   default='results/proteinmpnn',
                    help='Directory for job downloads (default: results/proteinmpnn)')
    ap.add_argument('--jobs-file',    default=None,
                    help='JSON file to save job IDs (default: jobs/<prefix>_mpnn_jobs.json)')
    args = ap.parse_args()

    pdb_path  = Path(args.pdb)
    out_root  = Path(args.output_dir) / args.prefix
    jobs_file = Path(args.jobs_file) if args.jobs_file else \
                Path('jobs') / f'{args.prefix}_mpnn_jobs.json'
    jobs_file.parent.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Find binding site ─────────────────────────────────────────────
    print(f"\n{'='*60}")
    if args.ligand:
        print(f" Finding {args.ligand} binding site in {pdb_path.name} chain {args.chain}")
        print(f"{'='*60}")
        iface_rnums = find_ligand_interface(pdb_path, args.chain, args.ligand, args.cutoff)
    elif args.partner_chain:
        print(f" Finding protein interface: chain {args.chain} vs chain {args.partner_chain}"
              f" in {pdb_path.name}")
        print(f"{'='*60}")
        iface_rnums = find_protein_interface(pdb_path, args.chain, args.partner_chain, args.cutoff)
    else:
        print(f" Using explicit interface residues for chain {args.chain} in {pdb_path.name}")
        print(f"{'='*60}")
        iface_rnums = sorted(int(x) for x in args.iface_residues.split(','))
    seq_positions = pdb_to_sequential(pdb_path, args.chain, iface_rnums)
    fixed_str     = f"{args.chain}:{','.join(map(str, seq_positions))}"
    print(f"  Interface residues ({len(iface_rnums)}) PDB nums: {iface_rnums}")
    print(f"  Sequential positions ({len(seq_positions)}):       {seq_positions}")
    print(f"  Fixed string: {fixed_str}")

    # ── Step 2: Submit jobs ───────────────────────────────────────────────────
    temps_outside  = np.round(np.linspace(args.temp_min, args.temp_max, args.n_outside),  2).tolist()
    temps_conserved= np.round(np.linspace(args.temp_min, args.temp_max, args.n_conserved), 2).tolist()

    print(f"\n Submitting {args.n_outside} mutoutside jobs (T={args.temp_min}–{args.temp_max})")
    print(f" Submitting {args.n_conserved} mutconserved jobs (T={args.temp_min}–{args.temp_max})")
    print()

    jobs = {}

    for i, t in enumerate(temps_outside, 1):
        name = f'{args.prefix}_mutoutside_{i}'
        jid  = submit_job(pdb_path, args.chain, fixed_str, t, name)
        if jid:
            jobs[name] = jid
            print(f"  [mutoutside  {i:>2}/{args.n_outside}] T={t:.2f}  job={jid[:8]}")

    for i, t in enumerate(temps_conserved, 1):
        name = f'{args.prefix}_mutconserved_{i}'
        jid  = submit_job(pdb_path, args.chain, None, t, name)
        if jid:
            jobs[name] = jid
            print(f"  [mutconserved {i:>2}/{args.n_conserved}] T={t:.2f}  job={jid[:8]}")

    with open(jobs_file, 'w') as f:
        json.dump(jobs, f, indent=2)
    print(f"\n  {len(jobs)} jobs submitted → {jobs_file}")

    # ── Step 3: Wait ──────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f" Waiting for {len(jobs)} jobs...")
    print(f"{'='*60}")
    wait_for_jobs(list(jobs.values()))

    # ── Step 4: Download and parse ────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f" Downloading and parsing results")
    print(f"{'='*60}")
    sequences = {}
    for name, jid in jobs.items():
        jdir = out_root / name
        download_job(jid, jdir)
        seq  = parse_fasta(jdir)
        if seq:
            sequences[name] = seq
            print(f"  {name}: {len(seq)} residues")
        else:
            print(f"  WARNING: no sequence parsed for {name}")

    # ── Step 5: Validation summary ────────────────────────────────────────────
    native_seqres = list(SeqIO.parse(str(pdb_path), 'pdb-seqres'))
    native_rec    = next((r for r in native_seqres
                          if r.id.endswith(f':{args.chain}')
                          or r.id.split(':')[-1] == args.chain), None)
    if native_rec:
        native_seq = str(native_rec.seq)
        fixed_idx  = [r - 1 for r in iface_rnums]
        print(f"\n{'Name':<40} {'T':>5} {'Total mut':>10} {'Site mut':>9} {'Conserved?':>11}")
        print('-' * 80)
        t_map = {f'{args.prefix}_mutoutside_{i+1}':  t for i,t in enumerate(temps_outside)}
        t_map.update({f'{args.prefix}_mutconserved_{i+1}': t for i,t in enumerate(temps_conserved)})
        for name, seq in sequences.items():
            t       = t_map.get(name, '?')
            tot_mut = sum(1 for a, b in zip(native_seq, seq) if a != b and b != 'X' and a != 'X')
            site_mut= sum(1 for i in fixed_idx if i < len(seq) and seq[i] != native_seq[i] and seq[i] != 'X')
            print(f"{name:<40} {t:>5} {tot_mut:>10} {site_mut:>9} {str(site_mut==0):>11}")

    # ── Step 6: Append to TSV ─────────────────────────────────────────────────
    total = append_to_tsv(args.output_tsv, sequences)
    print(f"\n  Added {len(sequences)} sequences → {args.output_tsv}  ({total} total rows)")


if __name__ == '__main__':
    main()
