"""
Fetch representative PDB structures from CATH superfamilies for structural_db.
Uses the CATH v4_3_0 cathtree API to list all domain representatives.
Files are saved as {label}_{domain_id}.pdb for clear identification.
"""

import os
import requests
import time

SUPERFAMILIES = {
    "3.90.1350.10": "translocation",
    "3.90.175.10": "catalytic",
    "2.60.120.200": "receptor_binding",
}

OUT_DIR = "structural_db"
os.makedirs(OUT_DIR, exist_ok=True)

CATH_API = "https://www.cathdb.info/version/v4_3_0/api/rest"


def get_domain_representatives(sfam_id):
    """Get representative domain IDs by traversing CATH tree to leaf depth (9)."""
    url = f"{CATH_API}/cathtree/from_cath_id_to_depth/{sfam_id}/9"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    data = r.json()

    domain_ids = []

    def traverse(node):
        if node.get("cath_id_depth") == 9:
            domain_id = node.get("example_domain_id")
            if domain_id:
                domain_ids.append(domain_id)
        for child in node.get("children", []):
            traverse(child)

    traverse(data)
    return domain_ids


def download_domain_pdb(domain_id, out_path):
    """Download domain coordinates from CATH."""
    url = f"{CATH_API}/id/{domain_id}.pdb"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        print(f"  WARNING: could not download {domain_id} (status {r.status_code})")
        return False
    with open(out_path, "w") as f:
        f.write(r.text)
    return True


for sfam_id, label in SUPERFAMILIES.items():
    print(f"\n=== {sfam_id} ({label}) ===")
    domain_dir = os.path.join(OUT_DIR, label)
    os.makedirs(domain_dir, exist_ok=True)

    try:
        reps = get_domain_representatives(sfam_id)
        print(f"  Found {len(reps)} domain representatives")
    except Exception as e:
        print(f"  ERROR fetching domain list: {e}")
        continue

    success = 0
    for domain_id in reps:
        # Label files as {label}_{domain_id}.pdb for clear identification
        out_path = os.path.join(domain_dir, f"{label}_{domain_id}.pdb")
        if os.path.exists(out_path):
            print(f"  {label}_{domain_id} already exists, skipping")
            success += 1
            continue
        print(f"  Downloading {domain_id}...")
        if download_domain_pdb(domain_id, out_path):
            success += 1
        time.sleep(0.2)  # be polite to the API

    print(f"  Done: {success}/{len(reps)} structures saved to {domain_dir}/")

print("\nAll done.")
