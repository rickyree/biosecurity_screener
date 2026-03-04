"""AlphaFold Database utility functions."""

import requests
import numpy as np
from pathlib import Path


def fetch_metadata(uniprot_id: str) -> dict:
    """Retrieve AlphaFold entry metadata for a UniProt accession."""
    endpoint = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    response = requests.get(endpoint)
    response.raise_for_status()
    return response.json()[0]


def download_alphafold_files(af_id: str, output_dir: str = "./") -> dict:
    """Download all files for an AlphaFold entry.

    Args:
        af_id: AlphaFold entry ID (e.g., 'AF-P04637-F1')
        output_dir: Directory to save files

    Returns:
        Dictionary with paths to downloaded files
    """
    base = "https://alphafold.ebi.ac.uk/files"
    files = {
        "structure": f"{af_id}-model_v4.cif",
        "confidence": f"{af_id}-confidence_v4.json",
        "pae": f"{af_id}-predicted_aligned_error_v4.json",
    }

    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    downloaded = {}
    for label, filename in files.items():
        url = f"{base}/{filename}"
        resp = requests.get(url)
        resp.raise_for_status()

        filepath = output_path / filename
        filepath.write_bytes(resp.content)
        downloaded[label] = str(filepath)

    return downloaded


def download_pdb(af_id: str, output_dir: str = "./") -> str:
    """Fetch structure in legacy PDB format."""
    url = f"https://alphafold.ebi.ac.uk/files/{af_id}-model_v4.pdb"
    resp = requests.get(url)
    resp.raise_for_status()

    outpath = Path(output_dir) / f"{af_id}.pdb"
    outpath.parent.mkdir(exist_ok=True)
    outpath.write_bytes(resp.content)
    return str(outpath)


def get_plddt_scores(af_id: str) -> dict:
    """Fetch and analyze per-residue confidence scores.

    Args:
        af_id: AlphaFold entry ID

    Returns:
        Dictionary with scores array and statistics
    """
    url = f"https://alphafold.ebi.ac.uk/files/{af_id}-confidence_v4.json"
    data = requests.get(url).json()
    scores = np.array(data['confidenceScore'])

    return {
        'scores': scores,
        'mean': float(np.mean(scores)),
        'median': float(np.median(scores)),
        'high_conf_count': int(np.sum(scores > 90)),
        'low_conf_count': int(np.sum(scores < 50)),
        'frac_high_conf': float(np.mean(scores > 90)),
        'frac_low_conf': float(np.mean(scores < 50)),
        'length': len(scores),
    }


def get_pae_matrix(af_id: str) -> np.ndarray:
    """Fetch predicted aligned error matrix.

    Args:
        af_id: AlphaFold entry ID

    Returns:
        2D numpy array of PAE values
    """
    url = f"https://alphafold.ebi.ac.uk/files/{af_id}-predicted_aligned_error_v4.json"
    data = requests.get(url).json()
    return np.array(data['distance'])


def plot_pae(af_id: str, output_path: str = None):
    """Generate heatmap of predicted aligned error.

    Args:
        af_id: AlphaFold entry ID
        output_path: Optional path to save figure

    Returns:
        matplotlib Figure object
    """
    import matplotlib.pyplot as plt

    pae = get_pae_matrix(af_id)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(pae, cmap='Greens_r', vmin=0, vmax=30)

    ax.set_xlabel('Scored Residue')
    ax.set_ylabel('Aligned Residue')
    ax.set_title(f'PAE Matrix: {af_id}')

    plt.colorbar(im, ax=ax, label='Expected Error (A)')

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')

    return fig


def compute_contact_map(cif_path: str, threshold: float = 8.0) -> dict:
    """Generate residue contact map from alpha carbon distances.

    Args:
        cif_path: Path to mmCIF file
        threshold: Distance cutoff in Angstroms

    Returns:
        Dictionary with distance matrix, contacts, and statistics
    """
    from Bio.PDB import MMCIFParser
    from scipy.spatial.distance import pdist, squareform

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_path)

    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_coords.append(residue['CA'].get_coord())

    coords = np.array(ca_coords)
    dist_matrix = squareform(pdist(coords))
    contacts = (dist_matrix > 0) & (dist_matrix < threshold)

    return {
        'distances': dist_matrix,
        'contacts': contacts,
        'num_residues': len(ca_coords),
        'num_contacts': int(np.sum(contacts) // 2),
    }


def batch_analyze(uniprot_ids: list) -> "pd.DataFrame":
    """Process multiple proteins and compile statistics.

    Args:
        uniprot_ids: List of UniProt accession IDs

    Returns:
        DataFrame with analysis results
    """
    import pandas as pd
    from Bio.PDB import alphafold_db

    records = []

    for uid in uniprot_ids:
        try:
            preds = list(alphafold_db.get_predictions(uid))
            if not preds:
                continue

            entry = preds[0]
            af_id = entry['entryId']
            stats = get_plddt_scores(af_id)

            records.append({
                'uniprot': uid,
                'alphafold_id': af_id,
                'length': stats['length'],
                'mean_plddt': stats['mean'],
                'median_plddt': stats['median'],
                'frac_high_conf': stats['frac_high_conf'],
                'frac_low_conf': stats['frac_low_conf'],
            })
        except Exception as e:
            print(f"Failed for {uid}: {e}")

    return pd.DataFrame(records)


def download_proteome(tax_id: int, output_dir: str = "./proteomes"):
    """Retrieve all AlphaFold predictions for a species by taxonomy ID.

    Args:
        tax_id: NCBI taxonomy ID (e.g., 9606 for human)
        output_dir: Directory to save files
    """
    import subprocess

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    gs_pattern = f"gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{tax_id}-*_v4.tar"

    subprocess.run(
        ["gsutil", "-m", "cp", gs_pattern, f"{output_dir}/"],
        check=True
    )
