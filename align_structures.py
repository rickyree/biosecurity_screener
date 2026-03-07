#!/usr/bin/env python
"""
Superimpose two PDB structures and calculate RMSD.
"""

import numpy as np
from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.PDB.vectors import Vector
import sys

def align_structures(pdb1_path, pdb2_path):
    """Align two structures and calculate RMSD."""
    parser = PDBParser(QUIET=True)

    # Load structures
    struct1 = parser.get_structure('struct1', pdb1_path)
    struct2 = parser.get_structure('struct2', pdb2_path)

    # Get first model and chain A from each
    model1 = struct1[0]
    model2 = struct2[0]

    # Try to find common chains
    chains1 = list(model1.get_chains())
    chains2 = list(model2.get_chains())

    print(f"Structure 1 chains: {[c.id for c in chains1]}")
    print(f"Structure 2 chains: {[c.id for c in chains2]}")

    # Use first chain from each (usually chain A)
    chain1 = chains1[0]
    chain2 = chains2[0]

    # Get CA atoms for alignment
    atoms1 = []
    atoms2 = []

    # Get residues with CA atoms
    residues1 = [r for r in chain1.get_residues() if r.id[0] == ' ' and 'CA' in r]
    residues2 = [r for r in chain2.get_residues() if r.id[0] == ' ' and 'CA' in r]

    print(f"Structure 1 residues: {len(residues1)}")
    print(f"Structure 2 residues: {len(residues2)}")

    # Align by sequence position (taking minimum length)
    min_len = min(len(residues1), len(residues2))

    for i in range(min_len):
        if 'CA' in residues1[i] and 'CA' in residues2[i]:
            atoms1.append(residues1[i]['CA'])
            atoms2.append(residues2[i]['CA'])

    print(f"Aligned {len(atoms1)} CA atoms")

    # Superimpose
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)

    rmsd = super_imposer.rms
    print(f"RMSD: {rmsd:.3f} Å")

    # Apply transformation
    super_imposer.apply(model2.get_atoms())

    return rmsd, super_imposer

if __name__ == '__main__':
    pdb1 = 'binding/2vsm.pdb'  # holo
    pdb2 = 'binding/2vwd_A_boltz2.pdb'  # apo

    print(f"Aligning {pdb1} chain A (holo) to {pdb2} (apo)")
    rmsd, superimposer = align_structures(pdb1, pdb2)

    print(f"\nFinal RMSD: {rmsd:.3f} Å")
    if rmsd < 2.0:
        print("Good alignment - structures are similar")
    elif rmsd < 5.0:
        print("Moderate alignment - some conformational differences")
    else:
        print("Poor alignment - significant structural differences")
