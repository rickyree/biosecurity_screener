# RCSB PDB API Reference

Technical reference for programmatic access to the Protein Data Bank.

## API Endpoints

| API | Purpose | Base URL |
|-----|---------|----------|
| Data API | Retrieve entry details by ID | `https://data.rcsb.org/rest/v1/` |
| Search API | Find entries matching criteria | `https://search.rcsb.org/rcsbsearch/v2/query` |
| GraphQL | Flexible structured queries | `https://data.rcsb.org/graphql` |
| File Server | Download coordinates | `https://files.rcsb.org/download/` |

## Data API Resources

| Resource | Endpoint | Example |
|----------|----------|---------|
| Entry | `/rest/v1/core/entry/{pdb_id}` | `4HHB` |
| Polymer Entity | `/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}` | `4HHB/1` |
| Nonpolymer Entity | `/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}` | `4HHB/5` |
| Assembly | `/rest/v1/core/assembly/{pdb_id}/{assembly_id}` | `4HHB/1` |
| Chemical Component | `/rest/v1/core/chem_comp/{comp_id}` | `HEM` |

### Common Data Fields

**Entry-level:**
- `struct.title` - Descriptive title
- `exptl[].method` - Experimental technique
- `rcsb_entry_info.resolution_combined` - Resolution (angstroms)
- `rcsb_entry_info.deposited_atom_count` - Atom count
- `rcsb_accession_info.deposit_date` - Submission date
- `rcsb_accession_info.initial_release_date` - Publication date

**Polymer entity:**
- `entity_poly.pdbx_seq_one_letter_code` - Amino acid sequence
- `rcsb_polymer_entity.formula_weight` - Molecular weight (Da)
- `rcsb_entity_source_organism.scientific_name` - Source species
- `rcsb_entity_source_organism.ncbi_taxonomy_id` - NCBI taxon ID

## Search Query Types

| Class | Description | Use Case |
|-------|-------------|----------|
| `TextQuery` | Full-text keyword search | General discovery |
| `AttributeQuery` | Filter by specific properties | Precise filtering |
| `SeqSimilarityQuery` | Sequence similarity (MMseqs2) | Homolog finding |
| `SeqMotifQuery` | Regex patterns in sequences | Motif search |
| `StructSimilarityQuery` | 3D fold similarity (BioZernike) | Structural homologs |
| `StructMotifQuery` | Geometric motif matching | Active site search |
| `ChemSimilarityQuery` | Ligand substructure search | Drug discovery |

## AttributeQuery Operators

| Operator | Type | Example |
|----------|------|---------|
| `exact_match` | String | `value="Homo sapiens"` |
| `contains_words` | String | `value="kinase inhibitor"` |
| `contains_phrase` | String | `value="ATP binding"` |
| `equals` | Numeric | `value=2.0` |
| `greater` | Numeric | `value=100` |
| `less` | Numeric | `value=2.5` |
| `range` | Numeric | `value=(1.5, 2.5)` |
| `exists` | Boolean | (field is populated) |
| `in` | List | `value=["X-RAY DIFFRACTION", "NEUTRON DIFFRACTION"]` |

## Searchable Attributes

Use string paths with `AttributeQuery`. Common attributes:

```python
# Resolution
"rcsb_entry_info.resolution_combined"

# Organism
"rcsb_entity_source_organism.scientific_name"
"rcsb_entity_source_organism.ncbi_taxonomy_id"

# Molecular weight
"rcsb_polymer_entity.formula_weight"

# Experimental method
"exptl.method"

# Quality metrics
"refine.ls_R_factor_R_free"

# Dates
"rcsb_accession_info.initial_release_date"
"rcsb_accession_info.deposit_date"

# Ligand codes
"rcsb_nonpolymer_entity_instance_container_identifiers.comp_id"
```

## Query Combination

```python
from rcsbapi.search import TextQuery, AttributeQuery

# AND: both conditions required
query = TextQuery("kinase") & AttributeQuery(...)

# OR: either condition
query = condition1 | condition2

# NOT: exclude matches
query = ~AttributeQuery(...)

# Complex example
query = (
    TextQuery("kinase") &
    AttributeQuery(
        attribute="rcsb_entity_source_organism.scientific_name",
        operator="exact_match",
        value="Homo sapiens"
    ) &
    ~AttributeQuery(
        attribute="rcsb_entry_info.resolution_combined",
        operator="greater",
        value=3.0
    )
)
```

## Sequence Search Parameters

```python
from rcsbapi.search import SeqSimilarityQuery

# Minimum 25 residues required
query = SeqSimilarityQuery(
    value="VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
    evalue_cutoff=1e-5,           # E-value threshold
    identity_cutoff=0.8,          # Minimum identity (0-1)
    sequence_type="protein",      # "protein", "dna", or "rna"
)
```

## Structure Search Options

```python
from rcsbapi.search import StructSimilarityQuery

# Search by entry
query = StructSimilarityQuery(
    structure_search_type="entry_id",
    entry_id="4HHB"
)

# Search by chain
query = StructSimilarityQuery(
    structure_search_type="chain_id",
    entry_id="4HHB",
    chain_id="A"
)

# Search by assembly
query = StructSimilarityQuery(
    structure_search_type="assembly_id",
    entry_id="4HHB",
    assembly_id="1"
)
```

## Controlling Output

```python
from rcsbapi.search import TextQuery

query = TextQuery("hemoglobin")

# PDB IDs (default)
ids = list(query())

# Limit results
ids = list(query(rows=100))

# Polymer entity IDs
entities = list(query(return_type="polymer_entity"))
# ['4HHB_1', '4HHB_2', ...]

# Assembly IDs
assemblies = list(query(return_type="assembly"))
# ['4HHB-1', '4HHB-2', ...]
```

## File Download URLs

| Content | URL Template |
|---------|--------------|
| PDB format | `https://files.rcsb.org/download/{ID}.pdb` |
| mmCIF format | `https://files.rcsb.org/download/{ID}.cif` |
| Structure factors | `https://files.rcsb.org/download/{ID}-sf.cif` |
| Assembly N | `https://files.rcsb.org/download/{ID}.pdb{N}` |
| FASTA | `https://www.rcsb.org/fasta/entry/{ID}` |

## Rate Limiting

Guidelines for responsible API usage:
- Start with 2-3 requests per second
- Implement exponential backoff on 429 responses
- Cache results to avoid redundant requests
- Use GraphQL to reduce request count for complex queries

```python
import time
import requests

def fetch_with_backoff(url: str, max_retries: int = 5) -> requests.Response:
    delay = 1.0
    for attempt in range(max_retries):
        resp = requests.get(url)
        if resp.status_code == 200:
            return resp
        elif resp.status_code == 429:
            time.sleep(delay)
            delay *= 2
        else:
            resp.raise_for_status()
    raise RuntimeError(f"Failed after {max_retries} attempts")
```

## HTTP Status Codes

| Code | Meaning | Action |
|------|---------|--------|
| 200 | Success | Process response |
| 404 | Entry not found | Check ID; may be obsoleted |
| 429 | Rate limited | Wait and retry with backoff |
| 500 | Server error | Retry after brief delay |

## Debugging

Inspect query JSON:
```python
from rcsbapi.search import TextQuery
import json

query = TextQuery("hemoglobin")
print(json.dumps(query.to_dict(), indent=2))
```

Test endpoints via curl:
```bash
# Verify entry exists
curl -I https://data.rcsb.org/rest/v1/core/entry/4HHB

# Test search API
curl 'https://search.rcsb.org/rcsbsearch/v2/query?json={"query":{"type":"terminal","service":"text","parameters":{"value":"test"}},"return_type":"entry"}'
```

## Additional Resources

- Data API Schema: https://data.rcsb.org/redoc/index.html
- GraphQL Explorer: https://data.rcsb.org/graphql
- Python Client Docs: https://rcsbapi.readthedocs.io/
- GitHub Issues: https://github.com/rcsb/py-rcsb-api/issues
