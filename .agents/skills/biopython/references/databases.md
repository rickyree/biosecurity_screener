# NCBI Database Access via Bio.Entrez

## Configuration

### Required Setup

```python
from Bio import Entrez

# MANDATORY: Identify yourself to NCBI
Entrez.email = "researcher@university.edu"
```

### Optional API Key

Register at https://www.ncbi.nlm.nih.gov/account/settings/ for increased rate limits:

```python
Entrez.api_key = "your_api_key_here"
```

| Configuration | Rate Limit |
|---------------|------------|
| Email only | 3 requests/second |
| Email + API key | 10 requests/second |

Biopython manages rate limiting automatically.

## Entrez Functions

### einfo: Database Information

```python
# List available databases
handle = Entrez.einfo()
result = Entrez.read(handle)
handle.close()
print(result["DbList"])

# Database statistics
handle = Entrez.einfo(db="protein")
info = Entrez.read(handle)
handle.close()
print(f"Records: {info['DbInfo']['Count']}")
```

### esearch: Search Databases

```python
# Basic search
handle = Entrez.esearch(db="nucleotide", term="BRCA1[Gene] AND human[Organism]")
results = Entrez.read(handle)
handle.close()

accessions = results["IdList"]
total = results["Count"]
```

**Search Parameters:**

| Parameter | Purpose | Example |
|-----------|---------|---------|
| `db` | Database name | `"pubmed"`, `"nucleotide"`, `"protein"` |
| `term` | Search query | `"insulin[Gene] AND mouse[Organism]"` |
| `retmax` | Max IDs returned | `100` |
| `retstart` | Starting index | `0` (for pagination) |
| `sort` | Sort order | `"relevance"`, `"pub_date"` |
| `reldate` | Days since publication | `365` (last year) |
| `datetype` | Date field | `"pdat"` (publication date) |

### esummary: Record Summaries

```python
handle = Entrez.esummary(db="pubmed", id="35758934,35690214")
summaries = Entrez.read(handle)
handle.close()

for doc in summaries:
    print(f"Title: {doc['Title']}")
    print(f"Authors: {doc['AuthorList']}")
```

### efetch: Retrieve Full Records

```python
# Fetch GenBank record
handle = Entrez.efetch(
    db="nucleotide",
    id="NM_007294",
    rettype="gb",
    retmode="text"
)
data = handle.read()
handle.close()

# Parse with SeqIO
from Bio import SeqIO
handle = Entrez.efetch(db="nucleotide", id="NM_007294", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
```

**Return Types by Database:**

| Database | rettype | Description |
|----------|---------|-------------|
| nucleotide | `"fasta"` | FASTA format |
| nucleotide | `"gb"` | GenBank flat file |
| protein | `"fasta"` | FASTA format |
| protein | `"gp"` | GenPept format |
| pubmed | `"abstract"` | Abstract text |
| pubmed | `"medline"` | MEDLINE format |

### elink: Cross-Database Links

```python
# Find proteins encoded by a nucleotide sequence
handle = Entrez.elink(dbfrom="nucleotide", db="protein", id="NM_007294")
links = Entrez.read(handle)
handle.close()

for linkset in links[0]["LinkSetDb"]:
    print(f"Link type: {linkset['LinkName']}")
    protein_ids = [link["Id"] for link in linkset["Link"]]
```

### epost: Upload ID Lists

Store large ID lists on NCBI servers for subsequent operations:

```python
ids = ["35758934", "35690214", "34567890"]
handle = Entrez.epost(db="pubmed", id=",".join(ids))
result = Entrez.read(handle)
handle.close()

query_key = result["QueryKey"]
web_env = result["WebEnv"]

# Use in later query
handle = Entrez.efetch(
    db="pubmed",
    query_key=query_key,
    WebEnv=web_env,
    rettype="medline"
)
```

### egquery: Search All Databases

```python
handle = Entrez.egquery(term="hemoglobin")
result = Entrez.read(handle)
handle.close()

for db_result in result["eGQueryResult"]:
    if int(db_result["Count"]) > 0:
        print(f"{db_result['DbName']}: {db_result['Count']} hits")
```

### espell: Spelling Suggestions

```python
handle = Entrez.espell(db="pubmed", term="bioinfromatcis")
correction = Entrez.read(handle)
handle.close()
print(f"Did you mean: {correction['CorrectedQuery']}?")
```

## Database-Specific Examples

### PubMed

```python
# Search for papers
handle = Entrez.esearch(db="pubmed", term="CRISPR gene editing", retmax=25)
results = Entrez.read(handle)
handle.close()

# Fetch abstracts
handle = Entrez.efetch(
    db="pubmed",
    id=results["IdList"],
    rettype="abstract",
    retmode="text"
)
abstracts = handle.read()
handle.close()
```

### GenBank / Nucleotide

```python
from Bio import SeqIO

# Search by gene and organism
handle = Entrez.esearch(
    db="nucleotide",
    term="cytochrome b[Gene] AND Panthera tigris[Organism]",
    retmax=5
)
results = Entrez.read(handle)
handle.close()

# Download sequences
for acc_id in results["IdList"]:
    handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print(f"{record.id}: {len(record)} bp")
```

### Protein

```python
# Search SwissProt entries
handle = Entrez.esearch(
    db="protein",
    term="hemoglobin[Protein Name] AND reviewed[Filter]"
)
results = Entrez.read(handle)
handle.close()

# Fetch GenPept format
handle = Entrez.efetch(
    db="protein",
    id=results["IdList"][:3],
    rettype="gp",
    retmode="text"
)
for record in SeqIO.parse(handle, "genbank"):
    print(f"{record.id}: {record.description}")
handle.close()
```

### Taxonomy

```python
# Look up organism
handle = Entrez.esearch(db="taxonomy", term="Escherichia coli")
results = Entrez.read(handle)
handle.close()

# Get taxonomic details
handle = Entrez.efetch(db="taxonomy", id=results["IdList"][0], retmode="xml")
taxon = Entrez.read(handle)
handle.close()

for record in taxon:
    print(f"Scientific name: {record['ScientificName']}")
    print(f"TaxID: {record['TaxId']}")
    print(f"Lineage: {record['Lineage']}")
```

### Gene

```python
# Search for gene records
handle = Entrez.esearch(db="gene", term="TP53[Gene] AND Homo sapiens[Organism]")
results = Entrez.read(handle)
handle.close()

# Fetch gene summary
handle = Entrez.efetch(db="gene", id=results["IdList"][0], retmode="xml")
gene_data = Entrez.read(handle)
handle.close()
```

## Batch Processing

### Paginated Downloads

```python
# Get total count
handle = Entrez.esearch(db="pubmed", term="machine learning", retmax=0)
result = Entrez.read(handle)
handle.close()
total = int(result["Count"])

# Download in batches
batch_size = 200
all_ids = []

for start in range(0, min(total, 1000), batch_size):
    handle = Entrez.esearch(
        db="pubmed",
        term="machine learning",
        retstart=start,
        retmax=batch_size
    )
    batch = Entrez.read(handle)
    handle.close()
    all_ids.extend(batch["IdList"])
```

### Using Search History

```python
# Search with history
handle = Entrez.esearch(
    db="nucleotide",
    term="mitochondrion[Title] AND complete genome",
    usehistory="y"
)
result = Entrez.read(handle)
handle.close()

web_env = result["WebEnv"]
query_key = result["QueryKey"]
count = int(result["Count"])

# Fetch using history (more efficient)
batch_size = 100
for start in range(0, count, batch_size):
    handle = Entrez.efetch(
        db="nucleotide",
        query_key=query_key,
        WebEnv=web_env,
        retstart=start,
        retmax=batch_size,
        rettype="fasta",
        retmode="text"
    )
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    # Process records
```

## Search Query Syntax

### Boolean Operators

```python
# AND - both terms required
term = "cancer AND immunotherapy"

# OR - either term
term = "HIV OR AIDS"

# NOT - exclude term
term = "virus NOT bacteriophage"
```

### Field Tags

| Tag | Field | Example |
|-----|-------|---------|
| `[Gene]` | Gene name | `"BRCA1[Gene]"` |
| `[Organism]` | Species | `"human[Organism]"` |
| `[Title]` | Title field | `"genome[Title]"` |
| `[Author]` | Author name | `"Smith J[Author]"` |
| `[Journal]` | Journal name | `"Nature[Journal]"` |
| `[PDAT]` | Publication date | `"2023[PDAT]"` |
| `[Filter]` | Special filters | `"reviewed[Filter]"` |

```python
# Complex query
term = '(BRCA1[Gene] OR BRCA2[Gene]) AND "Homo sapiens"[Organism] AND 2020:2024[PDAT]'
```

## Error Handling

```python
from urllib.error import HTTPError
import time

def safe_fetch(db, uid, rettype="fasta", max_retries=3):
    """Fetch with retry logic."""
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db=db, id=uid, rettype=rettype, retmode="text")
            data = handle.read()
            handle.close()
            return data
        except HTTPError as e:
            if e.code == 429:  # Rate limited
                time.sleep(2 ** attempt)
            elif e.code == 400:
                raise ValueError(f"Invalid ID: {uid}")
            else:
                raise
    raise RuntimeError(f"Failed after {max_retries} attempts")
```

## Recommendations

1. **Always set Entrez.email** before any operations
2. **Use API key** for intensive workflows (10x higher rate limit)
3. **Close handles** immediately after reading
4. **Cache downloaded data** locally to avoid redundant fetches
5. **Use search history** (WebEnv) for large result sets
6. **Batch operations** to minimize API calls
7. **Handle errors gracefully** - network issues are common
8. **Respect rate limits** - don't circumvent automatic throttling
9. **Validate IDs/accessions** before fetching
10. **Parse XML carefully** - structure varies between databases
