# bioRxiv REST API

Programmatic access to preprint metadata from the bioRxiv server.

## Base URL

```
https://api.biorxiv.org
```

## Endpoints

### Date Range Query

Retrieve preprints within a specified date window.

**Path:**
```
GET /details/biorxiv/{from_date}/{to_date}
GET /details/biorxiv/{from_date}/{to_date}/{category}
```

**Parameters:**

| Parameter | Format | Required | Description |
|-----------|--------|----------|-------------|
| `from_date` | YYYY-MM-DD | Yes | Query start date |
| `to_date` | YYYY-MM-DD | Yes | Query end date |
| `category` | string | No | Subject area filter |

**Sample Requests:**
```
https://api.biorxiv.org/details/biorxiv/2024-06-01/2024-06-30
https://api.biorxiv.org/details/biorxiv/2024-06-01/2024-06-30/bioinformatics
```

**Response:**
```json
{
  "messages": [
    {
      "status": "ok",
      "count": 87,
      "total": 87
    }
  ],
  "collection": [
    {
      "doi": "10.1101/2024.06.15.598234",
      "title": "Sample Preprint Title",
      "authors": "Kim J, Lee S, Park H",
      "author_corresponding": "Kim J",
      "author_corresponding_institution": "Seoul Institute of Technology",
      "date": "2024-06-15",
      "version": "1",
      "type": "new results",
      "license": "cc_by",
      "category": "bioinformatics",
      "jatsxml": "https://www.biorxiv.org/content/...",
      "abstract": "Abstract text...",
      "published": ""
    }
  ]
}
```

### DOI Lookup

Fetch metadata for a single preprint.

**Path:**
```
GET /details/biorxiv/{doi}
```

**Sample Request:**
```
https://api.biorxiv.org/details/biorxiv/10.1101/2024.06.15.598234
```

### Interval-Based Query

Retrieve recent preprints with pagination.

**Path:**
```
GET /pubs/biorxiv/{days_back}/{cursor}/{format}
```

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `days_back` | Number of days to look back |
| `cursor` | Pagination offset (start at 0, increment by 100) |
| `format` | Output format: `json` or `xml` |

**Sample Request:**
```
https://api.biorxiv.org/pubs/biorxiv/3/0/json
```

**Paginated Response:**
```json
{
  "messages": [
    {
      "status": "ok",
      "count": 100,
      "total": 342,
      "cursor": 100
    }
  ],
  "collection": [...]
}
```

## Record Schema

Each preprint in the `collection` array contains:

| Field | Type | Notes |
|-------|------|-------|
| `doi` | string | Unique identifier |
| `title` | string | Paper title |
| `authors` | string | Comma-separated list |
| `author_corresponding` | string | Contact author |
| `author_corresponding_institution` | string | Author affiliation |
| `date` | string | Post date (YYYY-MM-DD) |
| `version` | string | Revision number |
| `type` | string | "new results", "confirmatory", etc. |
| `license` | string | "cc_by", "cc_by_nc", etc. |
| `category` | string | Subject area |
| `abstract` | string | Summary text |
| `jatsxml` | string | Structured XML URL |
| `published` | string | Journal citation (blank if unpublished) |

## Subject Areas

Use these identifiers for category filtering:

**Life Sciences:**
- `animal-behavior-and-cognition`
- `biochemistry`
- `biophysics`
- `cell-biology`
- `developmental-biology`
- `ecology`
- `evolutionary-biology`
- `genetics`
- `genomics`
- `immunology`
- `microbiology`
- `molecular-biology`
- `neuroscience`
- `paleontology`
- `pathology`
- `pharmacology-and-toxicology`
- `physiology`
- `plant-biology`
- `zoology`

**Applied & Computational:**
- `bioengineering`
- `bioinformatics`
- `cancer-biology`
- `clinical-trials`
- `epidemiology`
- `synthetic-biology`
- `systems-biology`

**Other:**
- `scientific-communication-and-education`

## Full-Text Access

### PDF Retrieval

Construct PDF URLs from DOI and version:

```
https://www.biorxiv.org/content/{doi}v{version}.full.pdf
```

Example:
```
https://www.biorxiv.org/content/10.1101/2024.06.15.598234v1.full.pdf
```

### Web Page

```
https://www.biorxiv.org/content/{doi}v{version}
```

### JATS XML

Use the `jatsxml` field from API responses for structured markup.

## Usage Guidelines

**Rate Limits:**
- Minimum 500ms between requests
- Add delays for batch operations
- Implement exponential backoff on 429 responses

**Request Headers:**
```
User-Agent: YourApp/1.0 (contact@example.com)
```

**Result Caching:**
Save API responses locally to minimize redundant calls.

## HTTP Codes

| Status | Meaning |
|--------|---------|
| 200 | Success |
| 404 | Invalid DOI or endpoint |
| 429 | Rate limited - slow down |
| 500 | Server error - retry later |

## Code Examples

### Basic Request

```python
import requests
import time

API_BASE = "https://api.biorxiv.org"

def fetch_preprints(endpoint: str) -> dict:
    url = f"{API_BASE}/{endpoint}"
    headers = {"User-Agent": "ResearchBot/1.0"}

    resp = requests.get(url, headers=headers, timeout=30)
    resp.raise_for_status()

    time.sleep(0.5)
    return resp.json()

# Fetch immunology papers from July 2024
results = fetch_preprints("details/biorxiv/2024-07-01/2024-07-31/immunology")

for item in results["collection"][:5]:
    print(f"[{item['date']}] {item['title']}")
```

### Paginated Retrieval

```python
def fetch_recent_all(days: int) -> list:
    """Retrieve all preprints from the past N days."""
    papers = []
    cursor = 0

    while True:
        data = fetch_preprints(f"pubs/biorxiv/{days}/{cursor}/json")
        batch = data.get("collection", [])

        if not batch:
            break

        papers.extend(batch)

        meta = data["messages"][0]
        total = meta.get("total", 0)

        if cursor + len(batch) >= total:
            break

        cursor += 100

    return papers

# Get all preprints from last 7 days
week_papers = fetch_recent_all(7)
print(f"Retrieved {len(week_papers)} preprints")
```

### PDF Download

```python
def download_pdf(doi: str, version: str, save_path: str) -> bool:
    url = f"https://www.biorxiv.org/content/{doi}v{version}.full.pdf"

    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()

    with open(save_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=8192):
            f.write(chunk)

    return True

# Download a preprint PDF
download_pdf("10.1101/2024.06.15.598234", "1", "paper.pdf")
```

## External Links

- bioRxiv: https://www.biorxiv.org/
- API root: https://api.biorxiv.org/
- JATS spec: https://jats.nlm.nih.gov/
