#!/usr/bin/env python
"""
Compare sequences from PDB files.
"""

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def get_sequence_from_pdb(pdb_path, chain_id='A'):
    """Extract sequence from PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_path)

    # Get first model
    model = structure[0]

    # Get specified chain
    if chain_id not in model:
        print(f"Chain {chain_id} not found in {pdb_path}")
        return None

    chain = model[chain_id]

    # Get residues and convert to sequence
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    sequence = ''.join(seq1(r.resname) for r in residues)

    return sequence

def compare_sequences():
    # Compare the first pair
    print("=== Comparing first pair ===")
    seq1_a = get_sequence_from_pdb('binding/A_B_b2.pdb', 'A')
    seq1_b = get_sequence_from_pdb('binding/A_B_b3_A_unbound.pdb', 'A')

    print(f"A_B_b2.pdb chain A length: {len(seq1_a) if seq1_a else 'N/A'}")
    print(f"A_B_b3_A_unbound.pdb chain A length: {len(seq1_b) if seq1_b else 'N/A'}")

    if seq1_a and seq1_b:
        identical1 = seq1_a == seq1_b
        print(f"Sequences identical: {identical1}")
        if not identical1:
            # Count differences
            diffs = sum(1 for a, b in zip(seq1_a, seq1_b) if a != b)
            print(f"Number of differences: {diffs}/{min(len(seq1_a), len(seq1_b))}")

    # Compare the second pair
    print("\n=== Comparing second pair ===")
    seq2_a = get_sequence_from_pdb('binding/2vsm.pdb', 'A')
    seq2_b = get_sequence_from_pdb('binding/2vwd_A_boltz2.pdb', 'A')

    print(f"2vsm.pdb chain A length: {len(seq2_a) if seq2_a else 'N/A'}")
    print(f"2vwd_A_boltz2.pdb chain A length: {len(seq2_b) if seq2_b else 'N/A'}")

    if seq2_a and seq2_b:
        identical2 = seq2_a == seq2_b
        print(f"Sequences identical: {identical2}")
        if not identical2:
            # Count differences
            diffs = sum(1 for a, b in zip(seq2_a, seq2_b) if a != b)
            print(f"Number of differences: {diffs}/{min(len(seq2_a), len(seq2_b))}")

    # Cross-compare
    print("\n=== Cross comparison ===")
    if seq1_a and seq2_a:
        cross1 = seq1_a == seq2_a
        print(f"A_B_b2.pdb chain A == 2vsm.pdb chain A: {cross1}")

    if seq1_b and seq2_b:
        cross2 = seq1_b == seq2_b
        print(f"A_B_b3_A_unbound.pdb chain A == 2vwd_A_boltz2.pdb chain A: {cross2}")

    # Show first 50 characters for comparison if different
    if seq1_a and seq1_b and seq1_a != seq1_b:
        print(f"\nFirst 50 chars A_B_b2: {seq1_a[:50]}")
        print(f"First 50 chars A_B_b3: {seq1_b[:50]}")

    if seq2_a and seq2_b and seq2_a != seq2_b:
        print(f"\nFirst 50 chars 2vsm: {seq2_a[:50]}")
        print(f"First 50 chars 2vwd: {seq2_b[:50]}")

if __name__ == '__main__':
    compare_sequences()
