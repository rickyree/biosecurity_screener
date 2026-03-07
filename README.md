# Aegis: Structure-based Biosecurity Screening
A tool that screens user protein sequences to detect maliciously-redesigned ligand-binding toxins using interface signature detection and structural homology screening.

## Guided Workflow
### Installation
1. **Clone our repository**
   
   ```
   git clone https://github.com/rickyree/biosecurity_screener.git
   ```
   ```
   cd biosecurity_screener
   ```
2. **Create a conda environment**
   
    ```
    conda create -n aegis python=3.12
    conda activate aegis
    ```
    
3. **Install Amina CLI**

   Follow the installation instructions at: https://app.aminoanalytica.com/

### Running Aegis

4. **Run the Aegis screening pipeline**
   
   ```
   python main.py
   ```
   
### Structural homology validation

5. **Predict structures for flagged proteins**
   
    Use Amina CLI to run structure prediction (e.g., **Boltz-2**)

    Example:
   
   ```
   amina run boltz2 -s "XXX" -o results/boltz2
   ```
   

7. **Extend the structural database**

    - To add reference structures: Edit `scripts/get_cath_structures.py` and specify the desired **CATH superfamily ID**

    - Then download the structures into `structural_db`:
      ```
      python get_cath_structures.py
      ```
    
8. **Structural alignment**
   
    Use Amina CLI to calculate the **TM-scores** between the predicted structures and the reference structures in `structural_db`.

   Example:
   
   ```
   amina run usalign -m predicted.pdb -t reference.pdb -o results/usalign >> results/usalign/metrics.txt
   ```

9. **Obtain results**

   - Edit `scripts/summarise_usalign_metrics.py` and modify the correct input and output file paths
   - Edit `scripts/plot_tm_distributions.py` and modify the correct input and output file paths


    ```
    python summarise_usalign_metrics.py
    ```

    ```
    python plot_tm_distributions.py
    ```
   
