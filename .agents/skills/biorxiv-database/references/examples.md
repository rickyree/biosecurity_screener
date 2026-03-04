# Code Examples

Extended code examples for the bioRxiv Database skill.

## Table of Contents
- [Command-Line Examples](#command-line-examples)
- [Programmatic API Usage](#programmatic-api-usage)
- [Workflow Examples](#workflow-examples)
- [Data Pipeline Integration](#data-pipeline-integration)

## Command-Line Examples

### Keyword Searches

Basic keyword search with date range:

```bash
python scripts/biorxiv_client.py \
  --terms "protein folding" "AlphaFold" \
  --since 2024-03-01 \
  --until 2024-09-30 \
  --out folding_papers.json
```

Restrict to specific category:

```bash
python scripts/biorxiv_client.py \
  --terms "single cell" "RNA sequencing" \
  --recent 120 \
  --subject genomics \
  --out scrna_recent.json
```

Limit search scope to title only:

```bash
python scripts/biorxiv_client.py \
  --terms "transformer" \
  --fields title \
  --recent 180
```

### Author Lookups

Search by author name with date range:

```bash
python scripts/biorxiv_client.py \
  --author "Chen" \
  --since 2023-06-01 \
  --until 2024-06-01 \
  --out chen_publications.json
```

Default behavior (last 3 years when dates omitted):

```bash
python scripts/biorxiv_client.py \
  --author "Patel" \
  --out patel_recent.json
```

### Date-Based Retrieval

Fetch all preprints from a specific period:

```bash
python scripts/biorxiv_client.py \
  --since 2024-03-01 \
  --until 2024-03-31 \
  --out march_2024.json
```

Combine with category filtering:

```bash
python scripts/biorxiv_client.py \
  --since 2024-07-01 \
  --until 2024-07-31 \
  --subject immunology \
  --out immunology_july.json
```

Use relative time ranges:

```bash
python scripts/biorxiv_client.py \
  --recent 14 \
  --out past_two_weeks.json
```

### DOI-Based Lookup

Retrieve metadata for a specific preprint:

```bash
python scripts/biorxiv_client.py \
  --doi "10.1101/2024.05.22.594321" \
  --out paper_metadata.json
```

Full URLs are accepted:

```bash
python scripts/biorxiv_client.py \
  --doi "https://doi.org/10.1101/2024.05.22.594321"
```

### PDF Downloads

Download a preprint PDF:

```bash
python scripts/biorxiv_client.py \
  --doi "10.1101/2024.05.22.594321" \
  --fetch-pdf preprint.pdf
```

### Result Limits

```bash
python scripts/biorxiv_client.py \
  --terms "immunotherapy" \
  --recent 60 \
  --max 25 \
  --out immunotherapy_top25.json
```

## Programmatic API Usage

Import the client class for custom workflows:

```python
from scripts.biorxiv_client import PreprintClient

client = PreprintClient(debug=True)

# Keyword search
results = client.find_by_terms(
    terms=["enzyme engineering", "directed evolution"],
    since="2024-01-01",
    until="2024-12-31",
    subject="biochemistry"
)

# Author search
author_papers = client.find_by_author(
    name="Garcia",
    since="2023-01-01",
    until="2024-12-31"
)

# Direct DOI lookup
metadata = client.get_by_doi("10.1101/2024.05.22.594321")

# PDF retrieval
client.fetch_pdf(
    doi="10.1101/2024.05.22.594321",
    destination="enzyme_paper.pdf"
)

# Format output
formatted = client.normalize(metadata, include_abstract=True)
```

### Metadata-Only Retrieval

```python
from scripts.biorxiv_client import PreprintClient

client = PreprintClient()
papers = client.find_by_terms(terms=["proteomics"], since="2024-01-01", until="2024-12-31")
metadata_only = [client.normalize(p, include_abstract=False) for p in papers]
```

### Batch PDF Download

```python
import json
from scripts.biorxiv_client import PreprintClient

with open('search_results.json') as f:
    results = json.load(f)

client = PreprintClient(debug=True)

for idx, entry in enumerate(results['papers'][:5]):
    client.fetch_pdf(entry['doi'], f"downloads/paper_{idx+1}.pdf")
```

## Workflow Examples

### Systematic Review

```python
# Step 1: Broad search
# Run: python scripts/biorxiv_client.py \
#   --terms "drug delivery" "nanoparticles" \
#   --since 2022-01-01 \
#   --until 2024-12-31 \
#   --subject bioengineering \
#   --out nanoparticle_papers.json

# Step 2: Analyze results
import json

with open('nanoparticle_papers.json') as f:
    data = json.load(f)

print(f"Total matches: {data['count']}")

for paper in data['papers'][:3]:
    print(f"\n{paper['title']}")
    print(f"  Authors: {paper['authors']}")
    print(f"  Posted: {paper['posted']}")
    print(f"  DOI: {paper['doi']}")
```

```python
# Step 3: Download selected papers
from scripts.biorxiv_client import PreprintClient

client = PreprintClient()
target_dois = ["10.1101/2024.05.22.594321", "10.1101/2024.06.15.678901"]

for doi in target_dois:
    safe_name = doi.replace("/", "_").replace(".", "_") + ".pdf"
    client.fetch_pdf(doi, f"selected/{safe_name}")
```

### Publication Trend Analysis

```bash
python scripts/biorxiv_client.py \
  --terms "large language models" \
  --since 2021-01-01 \
  --until 2024-12-31 \
  --subject bioinformatics \
  --out llm_biology_trends.json
```

### Researcher Monitoring

```bash
# Monitor multiple researchers
for name in Garcia Thompson Lee; do
  python scripts/biorxiv_client.py \
    --author "$name" \
    --recent 180 \
    --out "${name}_preprints.json"
done
```

## Data Pipeline Integration

```python
import json
import pandas as pd

with open('search_results.json') as f:
    data = json.load(f)

df = pd.DataFrame(data['papers'])

print(f"Papers found: {len(df)}")
print(f"Date span: {df['posted'].min()} to {df['posted'].max()}")

print("\nMost frequent authors:")
author_counts = df['authors'].str.split(',').explode().str.strip().value_counts()
print(author_counts.head(5))

subset = df[df['posted'] >= '2024-06-01']
subset.to_csv('filtered_results.csv', index=False)
```

### Dynamic Date Ranges

```python
from datetime import datetime, timedelta

today = datetime.now()
quarter_ago = today - timedelta(days=91)

# Then run:
# python scripts/biorxiv_client.py \
#   --since {quarter_ago.strftime('%Y-%m-%d')} \
#   --until {today.strftime('%Y-%m-%d')}
```
