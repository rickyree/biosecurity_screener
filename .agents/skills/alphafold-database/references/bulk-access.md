# AlphaFold Bulk Data Access

## Google Cloud Storage

All AlphaFold predictions are available via GCS.

### Command Line (gsutil)

```bash
# Install gsutil
pip install gsutil

# Browse available data
gsutil ls gs://public-datasets-deepmind-alphafold-v4/

# Download human proteome (taxonomy ID 9606)
gsutil -m cp "gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-9606-*.tar" ./human/

# Get accession list
gsutil cp gs://public-datasets-deepmind-alphafold-v4/accession_ids.csv ./
```

### Common Taxonomy IDs

| Organism | Tax ID |
|----------|--------|
| Human | 9606 |
| Mouse | 10090 |
| E. coli | 83333 |
| Drosophila | 7227 |
| Yeast (S. cerevisiae) | 559292 |
| Arabidopsis | 3702 |

## BigQuery Metadata Queries

Query AlphaFold metadata at scale.

```python
from google.cloud import bigquery

def query_high_confidence_proteins(organism, min_confidence=0.8, limit=100):
    """Find proteins with high pLDDT fraction via BigQuery."""
    client = bigquery.Client()

    sql = f"""
    SELECT
        entryId,
        uniprotAccession,
        gene,
        globalMetricValue,
        fractionPlddtVeryHigh
    FROM `bigquery-public-data.deepmind_alphafold.metadata`
    WHERE organismScientificName = '{organism}'
        AND fractionPlddtVeryHigh > {min_confidence}
    ORDER BY fractionPlddtVeryHigh DESC
    LIMIT {limit}
    """

    return client.query(sql).to_dataframe()

# Example: Find high-quality human predictions
df = query_high_confidence_proteins("Homo sapiens", min_confidence=0.85)
```

### Available BigQuery Fields

- `entryId`: AlphaFold entry ID
- `uniprotAccession`: UniProt accession
- `gene`: Gene name
- `organismScientificName`: Species
- `globalMetricValue`: Overall pLDDT
- `fractionPlddtVeryHigh`: Fraction with pLDDT >90

## Batch Processing Example

```python
import pandas as pd
import numpy as np
from scripts.alphafold_utils import batch_analyze

# Analyze multiple proteins
uniprot_ids = ["P00519", "P04629", "P06239", "P12931"]
results = batch_analyze(uniprot_ids)
print(results.to_string(index=False))
```

## 3D-Beacons Federated API

Search across multiple structure databases.

```python
import requests

def query_3dbeacons(uniprot_id):
    """Search 3D-Beacons for available structure models."""
    url = f"https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/{uniprot_id}.json"
    resp = requests.get(url)
    data = resp.json()

    # Filter to AlphaFold entries only
    return [s for s in data.get('structures', []) if s['provider'] == 'AlphaFold DB']
```

## Links

- Google Cloud Console: https://console.cloud.google.com/marketplace/product/bigquery-public-data/deepmind-alphafold
- 3D-Beacons: https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/
