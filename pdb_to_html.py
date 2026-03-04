import py3Dmol
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import PDB
from Bio.SeqUtils import seq1
from anarcii import Anarcii


def save_pdb_as_html(pdb_path, output_html="molecule.html"):
    with open(pdb_path, 'r') as f:
        pdb_data = f.read()

    fmt = "mmcif" if pdb_path.endswith(('.cif', '.mmcif')) else "pdb"

    # Create the viewer object
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, fmt)
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()

    # Get the HTML string and save to a file
    html_content = view._make_html()
    with open(output_html, "w") as f:
        f.write(html_content)
    print(f"Interactive structure saved to {output_html}")


filepath = "/Users/hjlee/python_sandbox/biohack/results/esmfold/all_candidates/1f0l_A_mutoutside_1/esmfold_1f0l_A_mutoutside_1_structure.pdb" 
filename = os.path.basename(filepath)
save_pdb_as_html(filepath, filename.replace('.pdb', '.html'))