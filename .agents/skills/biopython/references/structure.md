# Protein Structure Analysis with Bio.PDB

## Structural Hierarchy (SMCRA Model)

Bio.PDB represents macromolecular structures as a hierarchy:

```
Structure
├── Model (0, 1, 2...)        # Multiple models for NMR structures
│   ├── Chain ("A", "B"...)   # Polypeptide chains
│   │   ├── Residue           # Amino acids, nucleotides, ligands
│   │   │   └── Atom          # Individual atoms
```

## Parsing Structure Files

### PDB Format

```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)  # Suppress warnings
structure = parser.get_structure("myprotein", "1abc.pdb")
```

### mmCIF Format

```python
from Bio.PDB import MMCIFParser

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("myprotein", "1abc.cif")
```

### Downloading from RCSB PDB

```python
from Bio.PDB import PDBList

pdb_list = PDBList()

# Download PDB format
pdb_list.retrieve_pdb_file("1ABC", file_format="pdb", pdir="./structures")

# Download mmCIF format (recommended)
pdb_list.retrieve_pdb_file("1ABC", file_format="mmCif", pdir="./structures")
```

## Navigating the Hierarchy

### Accessing Models

```python
# First model (most common)
model = structure[0]

# Iterate all models
for model in structure:
    print(f"Model {model.id}")
```

### Accessing Chains

```python
# Specific chain
chain_a = model["A"]

# All chains
for chain in model:
    print(f"Chain {chain.id}: {len(list(chain.get_residues()))} residues")
```

### Accessing Residues

```python
# All residues in chain
for residue in chain:
    print(f"{residue.resname} {residue.id[1]}")

# Specific residue (hetfield, sequence_number, insertion_code)
residue = chain[(" ", 42, " ")]  # Standard residue at position 42
```

**Residue ID tuple format:** `(hetfield, seqnum, icode)`
- `hetfield`: `" "` for standard, `"H_XXX"` for heteroatoms, `"W"` for water
- `seqnum`: Sequence number
- `icode`: Insertion code (usually `" "`)

### Accessing Atoms

```python
# All atoms in residue
for atom in residue:
    print(f"{atom.name}: {atom.coord}")

# Specific atom
ca_atom = residue["CA"]  # Alpha carbon
n_atom = residue["N"]    # Backbone nitrogen
```

### Shorthand Access

```python
# Direct chain access: structure[model][chain][residue][atom]
ca = structure[0]["A"][42]["CA"]

# Get all atoms in structure
all_atoms = list(structure.get_atoms())

# Get all residues
all_residues = list(structure.get_residues())
```

## Atom Properties

```python
atom = residue["CA"]

# Coordinates (NumPy array)
coord = atom.coord  # [x, y, z]

# Properties
element = atom.element        # "C"
bfactor = atom.bfactor        # Temperature factor
occupancy = atom.occupancy    # Occupancy (0-1)
serial = atom.serial_number   # PDB serial number
fullname = atom.fullname      # Full atom name with spacing
```

## Geometric Calculations

### Distance Between Atoms

```python
# Simple subtraction gives distance
distance = atom1 - atom2  # Returns float in Angstroms
```

### Angles

```python
from Bio.PDB.vectors import calc_angle

angle_rad = calc_angle(
    atom1.get_vector(),
    atom2.get_vector(),  # Vertex
    atom3.get_vector()
)
angle_deg = angle_rad * 180 / 3.14159
```

### Dihedral Angles

```python
from Bio.PDB.vectors import calc_dihedral

dihedral = calc_dihedral(
    atom1.get_vector(),
    atom2.get_vector(),
    atom3.get_vector(),
    atom4.get_vector()
)
```

## Secondary Structure (DSSP)

Requires DSSP executable installed on system.

```python
from Bio.PDB import DSSP

model = structure[0]
dssp = DSSP(model, "structure.pdb")

# Access results
for key in dssp:
    residue_id = key[1]
    ss_code = dssp[key][2]     # H, E, B, G, I, T, S, or -
    rel_asa = dssp[key][3]     # Relative solvent accessibility
    phi = dssp[key][4]         # Phi angle
    psi = dssp[key][5]         # Psi angle

    print(f"Res {residue_id}: {ss_code}, RSA={rel_asa:.2f}")
```

**Secondary structure codes:**

| Code | Structure |
|------|-----------|
| H | Alpha helix |
| E | Beta strand |
| B | Beta bridge |
| G | 3-10 helix |
| I | Pi helix |
| T | Turn |
| S | Bend |
| - | Coil/loop |

## Spatial Queries

### NeighborSearch for Proximity

```python
from Bio.PDB import NeighborSearch

atoms = list(structure.get_atoms())
ns = NeighborSearch(atoms)

# Atoms within radius
center = structure[0]["A"][100]["CA"].coord
nearby_atoms = ns.search(center, 5.0)  # 5 Angstrom radius

# Residues within radius
nearby_residues = ns.search(center, 5.0, level="R")

# Chains within radius
nearby_chains = ns.search(center, 10.0, level="C")
```

### Contact Map

```python
def compute_contact_map(chain, cutoff=8.0):
    """Generate residue contact map based on CA distances."""
    residues = [r for r in chain if r.has_id("CA")]
    n = len(residues)
    contacts = []

    for i in range(n):
        for j in range(i + 1, n):
            dist = residues[i]["CA"] - residues[j]["CA"]
            if dist < cutoff:
                contacts.append((i, j, dist))

    return contacts
```

## Structure Superimposition

```python
from Bio.PDB import Superimposer

# Get matching atoms from two structures
ref_atoms = [a for a in structure1[0]["A"].get_atoms() if a.name == "CA"]
mob_atoms = [a for a in structure2[0]["A"].get_atoms() if a.name == "CA"]

# Ensure equal length
min_len = min(len(ref_atoms), len(mob_atoms))
ref_atoms = ref_atoms[:min_len]
mob_atoms = mob_atoms[:min_len]

# Superimpose
sup = Superimposer()
sup.set_atoms(ref_atoms, mob_atoms)
sup.apply(structure2.get_atoms())

print(f"RMSD: {sup.rms:.2f} Angstroms")
```

## Writing Structure Files

### PDB Format

```python
from Bio.PDB import PDBIO

io = PDBIO()
io.set_structure(structure)
io.save("output.pdb")
```

### mmCIF Format

```python
from Bio.PDB import MMCIFIO

io = MMCIFIO()
io.set_structure(structure)
io.save("output.cif")
```

### Selective Output

```python
from Bio.PDB import Select

class ChainSelector(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

class BackboneSelector(Select):
    def accept_atom(self, atom):
        return atom.name in ["N", "CA", "C", "O"]

# Save only chain A
io = PDBIO()
io.set_structure(structure)
io.save("chain_a.pdb", ChainSelector(["A"]))

# Save only backbone
io.save("backbone.pdb", BackboneSelector())
```

## Sequence from Structure

```python
from Bio.PDB import Polypeptide

ppb = Polypeptide.PPBuilder()

for model in structure:
    for chain in model:
        for pp in ppb.build_peptides(chain):
            seq = pp.get_sequence()
            print(f"Chain {chain.id}: {seq}")
```

## Structure Analysis Functions

### Find Binding Site Residues

```python
def find_ligand_contacts(structure, ligand_chain, ligand_resnum, radius=4.5):
    """Find protein residues near a ligand."""
    from Bio.PDB import NeighborSearch

    # Get ligand atoms
    ligand = structure[0][ligand_chain][ligand_resnum]
    ligand_atoms = list(ligand.get_atoms())

    # Get protein atoms
    protein_atoms = []
    for chain in structure[0]:
        for residue in chain:
            if residue.id[0] == " ":  # Standard amino acid
                protein_atoms.extend(residue.get_atoms())

    # Search for contacts
    ns = NeighborSearch(protein_atoms)
    contact_residues = set()

    for lig_atom in ligand_atoms:
        nearby = ns.search(lig_atom.coord, radius, level="R")
        contact_residues.update(nearby)

    return list(contact_residues)
```

### Calculate Center of Mass

```python
import numpy as np

def center_of_mass(entity):
    """Calculate center of mass for structure entity."""
    atomic_masses = {"C": 12.0, "N": 14.0, "O": 16.0, "S": 32.0, "H": 1.0}

    coords = []
    masses = []

    for atom in entity.get_atoms():
        mass = atomic_masses.get(atom.element, 12.0)
        coords.append(atom.coord)
        masses.append(mass)

    coords = np.array(coords)
    masses = np.array(masses)

    return np.sum(coords * masses[:, np.newaxis], axis=0) / np.sum(masses)
```

### Ramachandran Data

```python
from Bio.PDB import Polypeptide

def extract_phi_psi(structure):
    """Extract backbone torsion angles."""
    angles = []

    for chain in structure[0]:
        polypeptides = Polypeptide.PPBuilder().build_peptides(chain)
        for pp in polypeptides:
            for residue, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
                if phi is not None and psi is not None:
                    angles.append({
                        'chain': chain.id,
                        'resnum': residue.id[1],
                        'resname': residue.resname,
                        'phi': phi,
                        'psi': psi
                    })

    return angles
```

### Check Missing Atoms

```python
def check_missing_backbone(structure):
    """Identify residues with missing backbone atoms."""
    backbone = ["N", "CA", "C", "O"]
    missing = []

    for residue in structure.get_residues():
        if residue.id[0] != " ":  # Skip heteroatoms
            continue

        for atom_name in backbone:
            if not residue.has_id(atom_name):
                missing.append({
                    'chain': residue.parent.id,
                    'resnum': residue.id[1],
                    'resname': residue.resname,
                    'missing': atom_name
                })

    return missing
```

## Recommendations

1. **Use mmCIF** for large structures and modern analyses
2. **Set QUIET=True** to suppress parser warnings in production
3. **Check for missing atoms** before geometric calculations
4. **Use NeighborSearch** for efficient spatial queries (not nested loops)
5. **Handle multiple models** appropriately for NMR structures
6. **Be aware of heteroatom residue IDs** (different format than standard)
7. **Validate structure quality** with DSSP or Ramachandran analysis
8. **Cache downloaded structures** locally
9. **Consider alternate conformations** - some residues have multiple positions
10. **Use Select classes** for targeted structure output
