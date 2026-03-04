# AlphaFold Confidence Score Interpretation

## pLDDT (Per-residue Local Distance Difference Test)

Stored in B-factor column of structure files.

| Score | Interpretation | Guidance |
|-------|---------------|----------|
| >90 | Very high | Reliable for atomic-level analysis |
| 70-90 | Confident | Backbone trustworthy, sidechains approximate |
| 50-70 | Low | Treat with caution; may indicate flexibility |
| <50 | Very low | Likely disordered or incorrectly predicted |

## PAE (Predicted Aligned Error)

Matrix showing inter-residue position confidence.

| Value | Meaning |
|-------|---------|
| <5 A | High confidence in relative domain positions |
| 5-10 A | Moderate certainty |
| >15 A | Uncertain arrangement; domains may be flexible |

## Practical Guidance

### When to Trust Predictions

- Use regions with pLDDT >70 for structural analysis
- Check PAE before assuming domain arrangements
- Low pLDDT often indicates intrinsic disorder (biologically relevant)

### When to Be Cautious

- Terminal regions frequently have low confidence
- Loop regions between domains
- Flexible linkers connecting structured domains

### Extracting pLDDT from Structure Files

```python
from Bio.PDB import MMCIFParser

def extract_bfactor_plddt(cif_path):
    """Read pLDDT values from B-factor column of mmCIF file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)

    residue_scores = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    residue_scores.append({
                        'residue': residue.get_resname(),
                        'position': residue.get_id()[1],
                        'plddt': residue['CA'].get_bfactor()
                    })
    return residue_scores
```

## References

- Jumper et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583-589.
- Varadi et al. (2024). AlphaFold Protein Structure Database in 2024. *Nucleic Acids Research* 52, D431-D438.
