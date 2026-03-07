# ESM2 Prefilter Standalone

Standalone ESM2-based sequence filtering for protein-protein interaction screening with integrated interface identification.

## Quick Usage

```bash
# Holo structure with ligand
python esm2_prefilter.py --structure complex_AB.pdb --bound-ligand LIG --candidates-tsv candidates.tsv

# Apo structure (automatically uses PeSTo)
python esm2_prefilter.py --structure protein_A.pdb --apo --candidates-tsv candidates.tsv

# Explicit residue list
python esm2_prefilter.py --structure protein_A.pdb --iface-residues "32,64,65,68" --candidates-tsv candidates.tsv
```

## Arguments

### Required
- `--structure PDB` - Input protein structure (indicate with --apo if apo)
- `--candidates-tsv TSV` - Candidates with `candidate` and `sequence_chain1` columns

### Interface Methods (auto-selected)
- **Default (holo)**: <5Å from ligand/partner in structure
- `--apo` - Structure is apo → automatically uses PeSTo predictions with Amina (highest confidence pocket)
- `--iface-residues "1,2,3"` - Explicit residue numbers (overrides above)

### Optional
- `--bound-chain ID` - Receptor chain in holo structure (auto-detect)
- `--bound-ligand RES` - Ligand residue name (for single-chain holo, otherwise will look for bound protein)
- `--out PATH` - Output TSV path (auto-generate)
- `--plddt-thresh FLOAT` - Filter by B-factor/pLDDT (default: 0)

## Input Format

**Candidates TSV**:
```
candidate	sequence_chain1
protein_1	MKLLVVVDEVHHGVGK...
protein_2	MPLVQGDKIKIFVGIK...
```

## Output

**Results TSV**: `results/esm2_prefilter_{structure_stem}[_apo|_explicit].tsv`

| Column | Description |
|--------|-------------|
| `candidate` | Sequence identifier |
| `seq_len` | Sequence length |
| `contact_score` | ESM2 inter-cluster contact score |
| `chem_sim` | Chemistry similarity score |

**Console Summary**:
```
Reference contact score: 0.0398
Contact threshold (50%): 0.0199
Chemistry sim threshold: 0.80
PASS: 180  |  FILTERED: 60
```

## Interface Methods

| Method | Trigger | Description |
|--------|---------|-------------|
| **5Å** | Holo structure (default) | Residues within 5Å of ligand/partner |
| **PeSTo** | `--apo` flag | ML predictions (top pocket only) |
| **Explicit** | `--iface-residues` | User-defined residue numbers |


## Notes

- ESM2 embeddings cached in `results/esm2_embed_cache/`
- PeSTo temp files in `results/interface_tmp/`
