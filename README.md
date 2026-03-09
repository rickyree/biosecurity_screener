# Aegis: Structure-based Biosecurity Screening
A tool that screens user protein sequences to detect maliciously-redesigned ligand-binding toxins using interface signature detection and structural homology screening (site: https://rickyree.github.io/slidedeck/).

## Guided Workflow
### Installation
1. **Clone our repository**
   
   ```
   git clone https://github.com/rickyree/biosecurity_screener.git
   ```
   ```
   cd biosecurity_screener
   ```
2. **Create a virtual environment and install requirements.txt**

    
3. **Install Amina CLI**

   Follow the installation instructions at: https://app.aminoanalytica.com/

### Running Aegis

4. **ESM2-based sequence filtering**
   
```bash
# Holo structure with ligand
python esm2_prefilter.py --structure binding/cholix.pdb --bound-ligand NAD --candidates-tsv candidate_sequences.tsv --out output.tsv

# Apo structure (automatically uses PeSTo)
python esm2_prefilter.py --structure binding/cholix_unbound.pdb --apo --candidates-tsv candidate_sequences.tsv --out output.tsv

# Explicit residue list
python esm2_prefilter.py --structure cholix.pdb --iface-residues "331, 335, 434, 435, 436, 437, 438, 441, 444, 445, 449, 450, 451, 452, 454, 455, 467, 468, 469, 470, 475, 476, 479, 483, 556, 557, 561" --candidates-tsv candidate_sequences.tsv --out output.tsv
```
#### Arguments

##### Required
- `--structure PDB` - Input protein structure (indicate with --apo if apo)
- `--candidates-tsv TSV` - Candidates with `candidate` and `sequence_chain1` columns

##### Interface Methods (auto-selected)
- **Default (holo)**: <5Å from ligand/partner in structure
- `--apo` - Structure is apo → automatically uses PeSTo predictions with Amina (highest confidence pocket)
- `--iface-residues "1,2,3"` - Explicit residue numbers (overrides above)

##### Optional
- `--bound-chain ID` - Receptor chain in holo structure (auto-detect)
- `--bound-ligand RES` - Ligand residue name (for single-chain holo, otherwise will look for bound protein)
- `--out PATH` - Output TSV path (default: ./output_{name of input PDB file}_{suffix}.tsv, suffix indicating if holo, apo or explicitly defined interface residues)
- `--plddt-thresh FLOAT` - Filter by B-factor/pLDDT (default: 0)

#### Input Format

**Candidates TSV**:
```
candidate	sequence_chain1
protein_1	MKLLVVVDEVHHGVGK...
protein_2	MPLVQGDKIKIFVGIK...
```

#### Output

**Results TSV**: stored as defined in --out PATH

| Column | Description |
|--------|-------------|
| `candidate` | Sequence identifier |
| `seq_len` | Sequence length |
| `chem_sim` | Aegis score |


#### Interface Methods

| Method | Trigger | Description |
|--------|---------|-------------|
| **5Å** | Holo structure (default) | Residues within 5Å of ligand/partner |
| **PeSTo** | `--apo` flag | ML predictions (top pocket only) |
| **Explicit** | `--iface-residues` | User-defined residue numbers |


#### Notes

- ESM2 embeddings cached in `results/esm2_embed_cache/`
- PeSTo temp files in `results/interface_tmp/`

   
### Structural homology validation

5. **Generate histogram**

    - Run to_run_histogram.py, with the filename variable specified as your output file.
    - **All will show as flagged by commec unless the candidate is present in to_run_results_bsspassed.tsv.**

7. **Predict structures for flagged proteins**
   
    Use Amina CLI to run structure prediction (e.g., **Boltz-2**)

    Example:
   
   ```
   amina run boltz2 -s "XXX" -o results/boltz2
   ```
   

8. **Extend the structural database**

    - To add reference structures: Edit `scripts/get_cath_structures.py` and specify the desired **CATH superfamily ID**

    - Then download the structures into `structural_db`:
      ```
      python scripts/fetch_cath_structures.py
      ```
    
9. **Structural alignment**
   
    Use Amina CLI to calculate the **TM-scores** between the predicted structures and the reference structures in `structural_db`.

   Example:
   
   ```
   amina run usalign -m predicted.pdb -t reference.pdb -o results/usalign >> results/usalign/metrics.txt
   ```

10. **Obtain results**

   Run the following scripts to summarise US-align metrics and generate visualisation of TM-score distribution.

    ```
    python scripts/summarise_usalign_metrics.py
    ```


    ```
    python scripts/plot_tm_distributions.py
    ```

## Repository Structure
   

```
biosecurity_screener/
│
├── scripts/                 # Helper scripts for downloading structures and analysing results
│   └── fetch_cath_structures.py
│   └── plot_tm_distributions.py
│   └── summarise_usalign_metrics.py
│
├── structural_db/           # Reference protein structures
│
├── results/                 # Output files and results
│
├── esm2_prefilter.py        # Main script to run Aegis
│
└── README.md
```
