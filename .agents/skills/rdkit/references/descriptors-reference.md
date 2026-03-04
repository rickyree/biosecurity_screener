# RDKit Molecular Descriptors

## Quick Start

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles('c1ccccc1O')

# Single descriptor
mw = Descriptors.MolWt(mol)

# All descriptors at once
props = Descriptors.CalcMolDescriptors(mol)
```

---

## Core Descriptors

### Mass Properties

| Descriptor | Function | Description |
|------------|----------|-------------|
| MolWt | `Descriptors.MolWt(mol)` | Average molecular weight |
| ExactMolWt | `Descriptors.ExactMolWt(mol)` | Monoisotopic mass |
| HeavyAtomMolWt | `Descriptors.HeavyAtomMolWt(mol)` | MW excluding hydrogens |

### Lipophilicity & Solubility

| Descriptor | Function | Description |
|------------|----------|-------------|
| MolLogP | `Descriptors.MolLogP(mol)` | Wildman-Crippen LogP |
| MolMR | `Descriptors.MolMR(mol)` | Molar refractivity |

### Surface Area

| Descriptor | Function | Description |
|------------|----------|-------------|
| TPSA | `Descriptors.TPSA(mol)` | Topological polar surface area |
| LabuteASA | `Descriptors.LabuteASA(mol)` | Labute approximate surface area |

### Hydrogen Bonding

| Descriptor | Function | Description |
|------------|----------|-------------|
| NumHDonors | `Descriptors.NumHDonors(mol)` | H-bond donors (N-H, O-H) |
| NumHAcceptors | `Descriptors.NumHAcceptors(mol)` | H-bond acceptors (N, O) |

### Atom Counts

| Descriptor | Function | Description |
|------------|----------|-------------|
| HeavyAtomCount | `Descriptors.HeavyAtomCount(mol)` | Non-hydrogen atoms |
| NumHeteroatoms | `Descriptors.NumHeteroatoms(mol)` | Non-C, non-H atoms |

### Ring Analysis

| Descriptor | Function | Description |
|------------|----------|-------------|
| RingCount | `Descriptors.RingCount(mol)` | Total ring count |
| NumAromaticRings | `Descriptors.NumAromaticRings(mol)` | Aromatic rings |
| NumSaturatedRings | `Descriptors.NumSaturatedRings(mol)` | Saturated rings |
| NumAliphaticRings | `Descriptors.NumAliphaticRings(mol)` | Non-aromatic rings |
| NumAromaticHeterocycles | `Descriptors.NumAromaticHeterocycles(mol)` | Aromatic rings with heteroatoms |

### Flexibility

| Descriptor | Function | Description |
|------------|----------|-------------|
| NumRotatableBonds | `Descriptors.NumRotatableBonds(mol)` | Rotatable single bonds |
| FractionCsp3 | `Descriptors.FractionCsp3(mol)` | Fraction of sp3 carbons |

### Drug-Likeness

| Descriptor | Function | Description |
|------------|----------|-------------|
| qed | `Descriptors.qed(mol)` | Quantitative Estimate of Drug-likeness |

---

## Complexity Measures

| Descriptor | Function | Description |
|------------|----------|-------------|
| BertzCT | `Descriptors.BertzCT(mol)` | Bertz complexity index |
| BalabanJ | `Descriptors.BalabanJ(mol)` | Branching descriptor |

---

## Shape Indices

### Kappa Shape

| Descriptor | Function |
|------------|----------|
| Kappa1 | `Descriptors.Kappa1(mol)` |
| Kappa2 | `Descriptors.Kappa2(mol)` |
| Kappa3 | `Descriptors.Kappa3(mol)` |

### Chi Connectivity

| Descriptor | Function |
|------------|----------|
| Chi0-Chi4 | `Descriptors.Chi0(mol)` through `Descriptors.Chi4(mol)` |
| Chi0v-Chi4v | `Descriptors.Chi0v(mol)` through `Descriptors.Chi4v(mol)` |

---

## Electronic Properties

### E-State Indices

| Descriptor | Function | Description |
|------------|----------|-------------|
| MaxEStateIndex | `Descriptors.MaxEStateIndex(mol)` | Maximum E-state |
| MinEStateIndex | `Descriptors.MinEStateIndex(mol)` | Minimum E-state |

### Partial Charges

| Descriptor | Function | Description |
|------------|----------|-------------|
| MaxPartialCharge | `Descriptors.MaxPartialCharge(mol)` | Most positive charge |
| MinPartialCharge | `Descriptors.MinPartialCharge(mol)` | Most negative charge |

---

## BCUT Descriptors

Burden-CAS-University of Texas eigenvalue descriptors:

| Descriptor | Function | Description |
|------------|----------|-------------|
| BCUT2D_MWHI | `Descriptors.BCUT2D_MWHI(mol)` | Highest eigenvalue (MW weighted) |
| BCUT2D_MWLOW | `Descriptors.BCUT2D_MWLOW(mol)` | Lowest eigenvalue (MW weighted) |
| BCUT2D_LOGPHI | `Descriptors.BCUT2D_LOGPHI(mol)` | Highest eigenvalue (LogP weighted) |
| BCUT2D_LOGPLOW | `Descriptors.BCUT2D_LOGPLOW(mol)` | Lowest eigenvalue (LogP weighted) |

---

## Common Applications

### Lipinski Rule of Five

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def check_lipinski(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    return {
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'passes': mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
    }
```

### Lead-Like Filter

```python
def is_leadlike(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotbonds = Descriptors.NumRotatableBonds(mol)

    return 250 <= mw <= 350 and logp <= 3.5 and rotbonds <= 7
```

### Drug-Likeness Profile

```python
def druglikeness_profile(mol):
    return {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotBonds': Descriptors.NumRotatableBonds(mol),
        'AromaticRings': Descriptors.NumAromaticRings(mol),
        'QED': Descriptors.qed(mol)
    }
```

---

## Batch Processing

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

# Calculate all descriptors for a list
def compute_descriptors(smiles_list):
    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            results.append(Descriptors.CalcMolDescriptors(mol))
        else:
            results.append(None)
    return results

# Get available descriptor names
names = [name for name, _ in Descriptors._descList]
```

---

## Descriptor Categories

| Category | Key Descriptors |
|----------|-----------------|
| Physicochemical | MolWt, MolLogP, MolMR, TPSA |
| Topological | BertzCT, BalabanJ, Kappa1-3 |
| Electronic | Partial charges, E-state indices |
| Connectivity | Chi indices (0-4, v, n variants) |
| Atom counts | HeavyAtomCount, NumHeteroatoms |
| Ring analysis | RingCount, NumAromaticRings |
| Drug-likeness | QED, Lipinski parameters |
| Flexibility | NumRotatableBonds, FractionCsp3 |

---

## Tips

1. **Batch calculation** with `CalcMolDescriptors()` is faster than individual calls
2. Some descriptors return `None` for invalid molecules
3. **Normalize** values for machine learning applications
4. **Select relevant** descriptors - not all 200+ are useful for every task
5. **3D descriptors** require conformer generation first
