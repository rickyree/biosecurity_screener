# RDKit API Reference

## rdkit.Chem - Core Module

### Molecule Creation

| Function | Description |
|----------|-------------|
| `MolFromSmiles(smiles, sanitize=True)` | Parse SMILES string |
| `MolFromSmarts(smarts)` | Parse SMARTS query pattern |
| `MolFromMolFile(path, sanitize=True, removeHs=True)` | Load from MOL file |
| `MolFromMolBlock(text, sanitize=True, removeHs=True)` | Parse MOL block string |
| `MolFromMol2File(path, sanitize=True, removeHs=True)` | Load from MOL2 file |
| `MolFromPDBFile(path, sanitize=True, removeHs=True)` | Load from PDB file |
| `MolFromInchi(inchi, sanitize=True, removeHs=True)` | Parse InChI identifier |

### Molecule Serialization

| Function | Description |
|----------|-------------|
| `MolToSmiles(mol, isomericSmiles=True, canonical=True)` | Generate SMILES |
| `MolToSmarts(mol, isomericSmarts=False)` | Generate SMARTS |
| `MolToMolBlock(mol, includeStereo=True, confId=-1)` | Generate MOL block |
| `MolToMolFile(mol, path, includeStereo=True, confId=-1)` | Save to MOL file |
| `MolToPDBBlock(mol, confId=-1)` | Generate PDB block |
| `MolToInchi(mol, options='')` | Generate InChI |
| `MolToInchiKey(mol, options='')` | Generate InChI key |

### Batch Processing

| Class/Function | Description |
|----------------|-------------|
| `SDMolSupplier(path, sanitize=True, removeHs=True)` | SDF file iterator |
| `ForwardSDMolSupplier(fileobj, sanitize=True)` | Forward-only SDF iterator (memory efficient) |
| `MultithreadedSDMolSupplier(path, numWriterThreads=1)` | Parallel SDF reader |
| `SmilesMolSupplier(path, delimiter=' ', titleLine=True)` | SMILES file iterator |
| `SDWriter(path)` | Write molecules to SDF |
| `SmilesWriter(path, delimiter=' ')` | Write SMILES file |

### Sanitization & Validation

| Function | Description |
|----------|-------------|
| `SanitizeMol(mol, sanitizeOps=SANITIZE_ALL)` | Validate and clean molecule |
| `DetectChemistryProblems(mol)` | Find validation issues |
| `AssignStereochemistry(mol, cleanIt=True)` | Assign stereo labels |
| `AssignStereochemistryFrom3D(mol, confId=-1)` | Infer stereo from 3D coords |

### Hydrogen Manipulation

| Function | Description |
|----------|-------------|
| `AddHs(mol, explicitOnly=False, addCoords=False)` | Add explicit hydrogens |
| `RemoveHs(mol, implicitOnly=False)` | Remove hydrogens |
| `RemoveAllHs(mol)` | Strip all hydrogens |

### Aromaticity

| Function | Description |
|----------|-------------|
| `SetAromaticity(mol, model=AROMATICITY_RDKIT)` | Apply aromaticity model |
| `Kekulize(mol, clearAromaticFlags=False)` | Convert to Kekule form |

### Fragment Operations

| Function | Description |
|----------|-------------|
| `GetMolFrags(mol, asMols=False)` | Get disconnected fragments |
| `FragmentOnBonds(mol, bondIndices, addDummies=True)` | Cut specific bonds |
| `ReplaceSubstructs(mol, query, replacement)` | Replace matching patterns |
| `DeleteSubstructs(mol, query)` | Remove matching patterns |

### Substructure Search

| Method | Description |
|--------|-------------|
| `mol.HasSubstructMatch(query)` | Check for match |
| `mol.GetSubstructMatch(query)` | Get first match indices |
| `mol.GetSubstructMatches(query, uniquify=True)` | Get all matches |

---

## Molecule Object Methods

### Atom Access

| Method | Description |
|--------|-------------|
| `mol.GetAtoms()` | Iterator over all atoms |
| `mol.GetAtomWithIdx(idx)` | Get atom by index |
| `mol.GetNumAtoms()` | Count atoms |
| `mol.GetNumHeavyAtoms()` | Count non-hydrogen atoms |

### Bond Access

| Method | Description |
|--------|-------------|
| `mol.GetBonds()` | Iterator over all bonds |
| `mol.GetBondWithIdx(idx)` | Get bond by index |
| `mol.GetNumBonds()` | Count bonds |

### Ring Information

| Method | Description |
|--------|-------------|
| `mol.GetRingInfo()` | Get ring info object |
| `GetSymmSSSR(mol)` | Smallest set of smallest rings |
| `ring_info.NumRings()` | Count rings |
| `ring_info.AtomRings()` | Atom indices per ring |

---

## Atom Object Methods

| Method | Description |
|--------|-------------|
| `atom.GetSymbol()` | Element symbol |
| `atom.GetAtomicNum()` | Atomic number |
| `atom.GetIdx()` | Atom index |
| `atom.GetDegree()` | Number of explicit bonds |
| `atom.GetFormalCharge()` | Formal charge |
| `atom.GetIsAromatic()` | Aromaticity flag |
| `atom.GetHybridization()` | SP, SP2, SP3, etc. |
| `atom.IsInRing()` | True if in any ring |
| `atom.IsInRingSize(n)` | True if in n-membered ring |
| `atom.GetChiralTag()` | Chirality specification |

---

## Bond Object Methods

| Method | Description |
|--------|-------------|
| `bond.GetBondType()` | SINGLE, DOUBLE, TRIPLE, AROMATIC |
| `bond.GetBeginAtomIdx()` | First atom index |
| `bond.GetEndAtomIdx()` | Second atom index |
| `bond.GetIsConjugated()` | Conjugation flag |
| `bond.GetIsAromatic()` | Aromaticity flag |
| `bond.IsInRing()` | True if in any ring |
| `bond.GetStereo()` | STEREONONE, STEREOZ, STEREOE |

---

## rdkit.Chem.AllChem - Extended Functions

### 2D Coordinates

| Function | Description |
|----------|-------------|
| `Compute2DCoords(mol)` | Generate 2D layout |
| `GenerateDepictionMatching2DStructure(mol, ref)` | Align to template |

### 3D Coordinates

| Function | Description |
|----------|-------------|
| `EmbedMolecule(mol, maxAttempts=0, randomSeed=-1)` | Generate single conformer |
| `EmbedMultipleConfs(mol, numConfs=10, randomSeed=-1)` | Generate multiple conformers |
| `ConstrainedEmbed(mol, core)` | Embed with constraints |

### Geometry Optimization

| Function | Description |
|----------|-------------|
| `UFFOptimizeMolecule(mol, maxIters=200, confId=-1)` | UFF minimization |
| `MMFFOptimizeMolecule(mol, maxIters=200, confId=-1)` | MMFF minimization |

### Conformer Analysis

| Function | Description |
|----------|-------------|
| `GetConformerRMS(mol, id1, id2)` | RMSD between conformers |
| `GetConformerRMSMatrix(mol)` | All-vs-all RMSD |
| `AlignMol(probe, ref)` | Align two molecules |
| `AlignMolConformers(mol)` | Align all conformers |

### Reactions

| Function | Description |
|----------|-------------|
| `ReactionFromSmarts(smarts)` | Create reaction |
| `reaction.RunReactants(reactants)` | Execute reaction |

---

## rdkit.Chem.rdFingerprintGenerator - Modern Fingerprint API

| Function | Description |
|----------|-------------|
| `GetMorganGenerator(radius=2, fpSize=2048)` | Morgan/circular fingerprints |
| `GetRDKitFPGenerator(minPath=1, maxPath=7, fpSize=2048)` | Topological fingerprints |
| `GetAtomPairGenerator(minDistance=1, maxDistance=30)` | Atom pair fingerprints |
| `GetTopologicalTorsionGenerator()` | Torsion fingerprints |
| `generator.GetFingerprint(mol)` | Generate bit vector |
| `generator.GetCountFingerprint(mol)` | Generate count vector |

---

## rdkit.Chem.MACCSkeys

| Function | Description |
|----------|-------------|
| `GenMACCSKeys(mol)` | Generate 166-bit MACCS keys |

---

## rdkit.DataStructs - Fingerprint Operations

### Similarity Metrics

| Function | Description |
|----------|-------------|
| `TanimotoSimilarity(fp1, fp2)` | Tanimoto coefficient |
| `DiceSimilarity(fp1, fp2)` | Dice coefficient |
| `CosineSimilarity(fp1, fp2)` | Cosine similarity |

### Bulk Operations

| Function | Description |
|----------|-------------|
| `BulkTanimotoSimilarity(fp, fps)` | Tanimoto vs list |
| `BulkDiceSimilarity(fp, fps)` | Dice vs list |

---

## rdkit.Chem.Draw - Visualization

| Function | Description |
|----------|-------------|
| `MolToImage(mol, size=(300,300), highlightAtoms=None)` | PIL image |
| `MolToFile(mol, path, size=(300,300))` | Save to file |
| `MolsToGridImage(mols, molsPerRow=3, subImgSize=(200,200))` | Grid image |
| `DrawMorganBit(mol, bitId, bitInfo)` | Show Morgan bit environment |

### Drawing Options (rdMolDraw2D.MolDrawOptions)

| Property | Description |
|----------|-------------|
| `.addAtomIndices` | Show atom indices |
| `.addStereoAnnotation` | Show stereochemistry |
| `.bondLineWidth` | Line thickness |

---

## rdkit.Chem.Scaffolds.MurckoScaffold

| Function | Description |
|----------|-------------|
| `GetScaffoldForMol(mol)` | Extract Murcko scaffold |
| `MakeScaffoldGeneric(mol)` | Generic scaffold (all C) |

---

## rdkit.Chem.rdMolHash - Molecular Hashing

| Hash Type | Description |
|-----------|-------------|
| `HashFunction.CanonicalSmiles` | Canonical SMILES |
| `HashFunction.MurckoScaffold` | Murcko scaffold hash |
| `HashFunction.Regioisomer` | Ignore stereochemistry |

```python
from rdkit.Chem import rdMolHash
hash_val = rdMolHash.MolHash(mol, rdMolHash.HashFunction.CanonicalSmiles)
```

---

## rdkit.Chem.MolStandardize - Structure Normalization

| Function | Description |
|----------|-------------|
| `Normalize(mol)` | Normalize functional groups |
| `Cleanup(mol)` | Full cleanup pipeline |
| `Uncharger().uncharge(mol)` | Neutralize charges |
| `TautomerEnumerator().Canonicalize(mol)` | Canonical tautomer |

---

## rdkit.ML.Cluster.Butina - Clustering

| Function | Description |
|----------|-------------|
| `ClusterData(distances, nPts, distThresh, isDistData=True)` | Butina clustering |

---

## Constants

### Sanitization Flags

| Flag | Description |
|------|-------------|
| `SANITIZE_NONE` | Skip all |
| `SANITIZE_ALL` | All operations (default) |
| `SANITIZE_KEKULIZE` | Kekulize aromatics |
| `SANITIZE_SETAROMATICITY` | Perceive aromaticity |

### Bond Types

| Type | Description |
|------|-------------|
| `BondType.SINGLE` | Single bond |
| `BondType.DOUBLE` | Double bond |
| `BondType.TRIPLE` | Triple bond |
| `BondType.AROMATIC` | Aromatic bond |

### Hybridization

| Type | Description |
|------|-------------|
| `HybridizationType.SP` | SP hybridization |
| `HybridizationType.SP2` | SP2 hybridization |
| `HybridizationType.SP3` | SP3 hybridization |

### Chirality

| Type | Description |
|------|-------------|
| `ChiralType.CHI_UNSPECIFIED` | Not specified |
| `ChiralType.CHI_TETRAHEDRAL_CW` | Clockwise |
| `ChiralType.CHI_TETRAHEDRAL_CCW` | Counter-clockwise |
