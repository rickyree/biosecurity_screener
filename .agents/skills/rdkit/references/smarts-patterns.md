# SMARTS Pattern Reference

## Quick Start

```python
from rdkit import Chem

# Create pattern
pattern = Chem.MolFromSmarts('[CH2][OH1]')

# Search molecule
mol = Chem.MolFromSmiles('CCO')
matches = mol.GetSubstructMatches(pattern)
```

---

## Functional Groups

### Alcohols

| Type | SMARTS |
|------|--------|
| Primary alcohol | `[CH2][OH1]` |
| Secondary alcohol | `[CH1]([OH1])[C]` |
| Tertiary alcohol | `[C]([OH1])([C])([C])[C]` |
| Any alcohol | `[OH1][C]` |
| Phenol | `c[OH1]` |

### Carbonyl Compounds

| Type | SMARTS |
|------|--------|
| Aldehyde | `[CH1](=O)` |
| Ketone | `[C](=O)[C]` |
| Any carbonyl | `[C](=O)` |

### Carboxylic Acids & Derivatives

| Type | SMARTS |
|------|--------|
| Carboxylic acid | `C(=O)[OH1]` |
| Ester | `C(=O)O[C]` |
| Amide | `C(=O)N` |
| Acyl chloride | `C(=O)Cl` |
| Anhydride | `C(=O)OC(=O)` |

### Amines

| Type | SMARTS |
|------|--------|
| Primary amine | `[NH2][C]` |
| Secondary amine | `[NH1]([C])[C]` |
| Tertiary amine | `[N]([C])([C])[C]` |
| Aniline | `c[NH2]` |
| Any amine | `[NX3]` |

### Ethers & Halides

| Type | SMARTS |
|------|--------|
| Aliphatic ether | `[C][O][C]` |
| Any alkyl halide | `[C][F,Cl,Br,I]` |
| Aryl halide | `c[F,Cl,Br,I]` |

### Nitrogen & Sulfur Groups

| Type | SMARTS |
|------|--------|
| Nitrile | `C#N` |
| Nitro group | `[N+](=O)[O-]` |
| Thiol | `[C][SH1]` |
| Sulfide | `[C][S][C]` |
| Disulfide | `[C][S][S][C]` |
| Sulfone | `[C][S](=O)(=O)[C]` |

---

## Ring Systems

### Simple Rings

| Type | SMARTS |
|------|--------|
| Benzene | `c1ccccc1` |
| Cyclohexane | `C1CCCCC1` |
| 5-membered ring | `[r5]` |
| 6-membered ring | `[r6]` |
| Any ring atom | `[R]` |

### Heterocycles

| Type | SMARTS |
|------|--------|
| Pyridine | `n1ccccc1` |
| Pyrrole | `n1cccc1` |
| Furan | `o1cccc1` |
| Thiophene | `s1cccc1` |
| Imidazole | `n1cncc1` |

### Fused Rings

| Type | SMARTS |
|------|--------|
| Naphthalene | `c1ccc2ccccc2c1` |
| Indole | `c1ccc2[nH]ccc2c1` |
| Quinoline | `n1cccc2ccccc12` |

---

## Pharmacophore Features

### Hydrogen Bond Donors

| Type | SMARTS |
|------|--------|
| Hydroxyl | `[OH]` |
| Amine | `[NH,NH2]` |
| Any donor | `[OH,NH,NH2,NH3+]` |

### Hydrogen Bond Acceptors

| Type | SMARTS |
|------|--------|
| Carbonyl oxygen | `[O]=[C,S,P]` |
| Ether oxygen | `[OX2]` |
| Any acceptor | `[O,N]` |

---

## Drug-Like Scaffolds

| Type | SMARTS |
|------|--------|
| Biphenyl | `c1ccccc1-c2ccccc2` |
| Piperazine | `N1CCNCC1` |
| Piperidine | `N1CCCCC1` |
| Morpholine | `N1CCOCC1` |
| Benzamide | `c1ccccc1C(=O)N` |
| Sulfonamide | `S(=O)(=O)N` |
| Urea | `[N][C](=O)[N]` |

---

## Reactive & Toxicity Groups

### Electrophiles

| Type | SMARTS |
|------|--------|
| Acyl chloride | `C(=O)Cl` |
| Epoxide | `C1OC1` |
| Michael acceptor | `C=C[C](=O)` |

### Toxicity Alerts

| Type | SMARTS |
|------|--------|
| Catechol | `c1ccc(O)c(O)c1` |
| Quinone | `O=C1C=CC(=O)C=C1` |
| Reactive halide | `[C][I,Br]` |

---

## Atom & Bond Specifications

### Atom Types

| Type | SMARTS |
|------|--------|
| Any atom | `[*]` |
| Carbon (any) | `[C,c]` |
| Aliphatic carbon | `[C]` |
| Aromatic carbon | `[c]` |
| Heteroatom | `[!C;!H]` |
| Halogen | `[F,Cl,Br,I]` |

### Connectivity

| Type | SMARTS |
|------|--------|
| Terminal atom | `[D1]` |
| 2 neighbors | `[D2]` |
| Branch point | `[D3]` |
| Methyl group | `[CH3]` |
| SP3 carbon | `[CX4]` |
| SP2 carbon | `[CX3]` |

### Ring Membership

| Type | SMARTS |
|------|--------|
| In any ring | `[R]` |
| Not in ring | `[!R]` |
| In exactly 1 ring | `[R1]` |
| In 5-7 membered ring | `[r{5-7}]` |

### Charge & Stereo

| Type | SMARTS |
|------|--------|
| Positive charge | `[+]` |
| Negative charge | `[-]` |
| Clockwise chiral | `[C@]` |
| Counterclockwise | `[C@@]` |

---

## Usage Examples

### Check for Pattern

```python
from rdkit import Chem

patterns = {
    'alcohol': '[OH1][C]',
    'amine': '[NH2,NH1][C]',
    'carboxylic_acid': 'C(=O)[OH1]'
}

mol = Chem.MolFromSmiles('OCC(=O)O')

for name, smarts in patterns.items():
    query = Chem.MolFromSmarts(smarts)
    if mol.HasSubstructMatch(query):
        print(f"Contains {name}")
```

### Count Matches

```python
pattern = Chem.MolFromSmarts('[OH1]')
matches = mol.GetSubstructMatches(pattern)
print(f"Found {len(matches)} hydroxyl groups")
```

### Filter Compounds

```python
def filter_by_pattern(smiles_list, smarts):
    query = Chem.MolFromSmarts(smarts)
    return [
        smi for smi in smiles_list
        if Chem.MolFromSmiles(smi) and
           Chem.MolFromSmiles(smi).HasSubstructMatch(query)
    ]
```

---

## SMARTS Syntax Summary

| Symbol | Meaning |
|--------|---------|
| `[C]` | Aliphatic carbon |
| `[c]` | Aromatic carbon |
| `[CX4]` | Carbon with 4 connections |
| `[CH3]` | Methyl group |
| `[R]` | In ring |
| `[r6]` | In 6-membered ring |
| `[D2]` | Degree 2 (2 neighbors) |
| `[+]` | Positive charge |
| `[!C]` | Not carbon |
| `[#6]` | Atomic number 6 (carbon) |
| `~` | Any bond type |
| `-` | Single bond |
| `=` | Double bond |
| `#` | Triple bond |
| `:` | Aromatic bond |

---

## Tips

1. **Lowercase = aromatic** - `c`, `n`, `o`, `s` are aromatic atoms
2. **Use brackets for clarity** - `[C]` vs `C` have different meanings
3. **Start specific, then generalize** - begin with strict patterns
4. **Test patterns** - validate on known molecules first
5. **Check ring membership** - use `[R]` and `[r{n}]` appropriately
6. **Consider hydrogens** - implicit unless specified with `[H]` or counts
