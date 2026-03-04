# BLAST Operations

## Remote BLAST via NCBI

### Running Searches

```python
from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Load query sequence
query = SeqIO.read("query.fasta", "fasta")

# Execute BLAST
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="swissprot",
    sequence=str(query.seq)
)

# Save results
with open("results.xml", "w") as output:
    output.write(result_handle.read())
result_handle.close()
```

### BLAST Programs

| Program | Query | Database | Use Case |
|---------|-------|----------|----------|
| `blastn` | Nucleotide | Nucleotide | DNA/RNA similarity |
| `blastp` | Protein | Protein | Protein homology |
| `blastx` | Nucleotide | Protein | Find protein hits for DNA |
| `tblastn` | Protein | Nucleotide | Search DNA with protein |
| `tblastx` | Nucleotide | Nucleotide | 6-frame translation comparison |

### Common Databases

**Nucleotide:**
- `nt` - Non-redundant nucleotide (comprehensive)
- `refseq_rna` - RefSeq RNA sequences

**Protein:**
- `nr` - Non-redundant protein (comprehensive)
- `swissprot` - Curated UniProtKB/Swiss-Prot
- `refseq_protein` - RefSeq proteins
- `pdb` - Protein Data Bank sequences

### Search Parameters

```python
result_handle = NCBIWWW.qblast(
    program="blastn",
    database="nt",
    sequence=query_seq,
    expect=0.01,           # E-value threshold
    hitlist_size=100,      # Max hits
    alignments=50,         # Max alignments shown
    word_size=11,          # Initial word match size
    gapcosts="5 2",        # Gap open/extend penalties
    format_type="XML",     # Output format
    entrez_query="Homo sapiens[Organism]"  # Restrict to organism
)
```

### Query Input Options

```python
# Sequence string
NCBIWWW.qblast("blastn", "nt", "ATCGATCGATCGATCG")

# FASTA string
with open("query.fasta") as f:
    NCBIWWW.qblast("blastn", "nt", f.read())

# GenBank accession
NCBIWWW.qblast("blastn", "nt", "NM_007294")
```

## Parsing BLAST Results

### XML Output (Recommended)

```python
from Bio.Blast import NCBIXML

# Single result
with open("results.xml") as handle:
    record = NCBIXML.read(handle)

# Multiple results
with open("batch_results.xml") as handle:
    for record in NCBIXML.parse(handle):
        process_record(record)
```

### Accessing Result Data

```python
# Query information
print(f"Query: {record.query}")
print(f"Query length: {record.query_length}")
print(f"Database: {record.database}")
print(f"Database size: {record.database_sequences}")

# Iterate through hits
for alignment in record.alignments:
    print(f"\nHit: {alignment.title}")
    print(f"Accession: {alignment.accession}")
    print(f"Length: {alignment.length}")

    # High-scoring pairs (HSPs)
    for hsp in alignment.hsps:
        print(f"  E-value: {hsp.expect}")
        print(f"  Bit score: {hsp.bits}")
        print(f"  Identity: {hsp.identities}/{hsp.align_length}")
        print(f"  Gaps: {hsp.gaps}")
        print(f"  Query range: {hsp.query_start}-{hsp.query_end}")
        print(f"  Subject range: {hsp.sbjct_start}-{hsp.sbjct_end}")
```

### HSP Attributes

| Attribute | Description |
|-----------|-------------|
| `expect` | E-value |
| `score` | Raw alignment score |
| `bits` | Bit score |
| `identities` | Number of identical positions |
| `positives` | Similar positions (proteins) |
| `gaps` | Gap count |
| `align_length` | Alignment length |
| `query` | Aligned query sequence |
| `sbjct` | Aligned subject sequence |
| `match` | Match line showing identities |
| `query_start`, `query_end` | Query coordinates |
| `sbjct_start`, `sbjct_end` | Subject coordinates |

### Filtering Results

```python
E_THRESHOLD = 1e-5
IDENTITY_THRESHOLD = 0.7

for alignment in record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect > E_THRESHOLD:
            continue

        identity = hsp.identities / hsp.align_length
        if identity < IDENTITY_THRESHOLD:
            continue

        print(f"{alignment.accession}: {identity:.1%} identity, E={hsp.expect:.2e}")
```

## Local BLAST

### Prerequisites

Install BLAST+ command-line tools from NCBI.

### Creating Databases

```python
from Bio.Blast.Applications import NcbimakeblastdbCommandline

make_db = NcbimakeblastdbCommandline(
    input_file="sequences.fasta",
    dbtype="nucl",  # or "prot"
    out="my_database"
)
stdout, stderr = make_db()
```

### Running Local Searches

```python
from Bio.Blast.Applications import NcbiblastnCommandline

blastn = NcbiblastnCommandline(
    query="query.fasta",
    db="my_database",
    evalue=0.001,
    outfmt=5,  # XML format
    out="local_results.xml"
)
stdout, stderr = blastn()

# Parse results
from Bio.Blast import NCBIXML
with open("local_results.xml") as handle:
    record = NCBIXML.read(handle)
```

### Command-Line Wrappers

| Wrapper | Program |
|---------|---------|
| `NcbiblastnCommandline` | blastn |
| `NcbiblastpCommandline` | blastp |
| `NcbiblastxCommandline` | blastx |
| `NcbitblastnCommandline` | tblastn |
| `NcbitblastxCommandline` | tblastx |

## Tabular Output

### Format 6 (Tab-delimited)

```python
from Bio.Blast.Applications import NcbiblastnCommandline

blastn = NcbiblastnCommandline(
    query="query.fasta",
    db="database",
    outfmt=6,
    out="results.tsv"
)
stdout, stderr = blastn()

# Parse tabular output
with open("results.tsv") as f:
    for line in f:
        fields = line.strip().split('\t')
        query_id = fields[0]
        subject_id = fields[1]
        pct_identity = float(fields[2])
        alignment_length = int(fields[3])
        evalue = float(fields[10])
        bitscore = float(fields[11])
```

**Default columns (outfmt 6):**
qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore

### Custom Column Selection

```python
blastn = NcbiblastnCommandline(
    query="query.fasta",
    db="database",
    outfmt="6 qseqid sseqid pident evalue bitscore qseq sseq",
    out="custom_results.tsv"
)
```

## Analysis Utilities

### Extract Best Hits

```python
def get_top_hits(blast_record, n=10, max_evalue=0.001):
    """Extract top N hits below E-value threshold."""
    hits = []
    for alignment in blast_record.alignments[:n]:
        if alignment.hsps:
            best_hsp = alignment.hsps[0]
            if best_hsp.expect <= max_evalue:
                hits.append({
                    'accession': alignment.accession,
                    'title': alignment.title,
                    'evalue': best_hsp.expect,
                    'identity': best_hsp.identities / best_hsp.align_length,
                    'coverage': best_hsp.align_length / blast_record.query_length
                })
    return hits
```

### Calculate Percent Identity

```python
def percent_identity(hsp):
    """Calculate percent identity from HSP."""
    return (hsp.identities / hsp.align_length) * 100
```

### Fetch Hit Sequences

```python
from Bio import Entrez, SeqIO

def fetch_hit_sequences(blast_record, n=5):
    """Download sequences for top BLAST hits."""
    Entrez.email = "researcher@institution.edu"

    sequences = []
    for alignment in blast_record.alignments[:n]:
        handle = Entrez.efetch(
            db="protein",
            id=alignment.accession,
            rettype="fasta",
            retmode="text"
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequences.append(record)

    return sequences
```

## Batch BLAST

### Multiple Queries

```python
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

queries = list(SeqIO.parse("queries.fasta", "fasta"))

results = []
for query in queries:
    print(f"BLASTing {query.id}...")
    handle = NCBIWWW.qblast("blastp", "swissprot", str(query.seq))
    record = NCBIXML.read(handle)
    handle.close()
    results.append((query.id, record))

# Summarize results
for query_id, record in results:
    n_hits = len(record.alignments)
    best_e = record.alignments[0].hsps[0].expect if record.alignments else None
    print(f"{query_id}: {n_hits} hits, best E={best_e}")
```

## Reciprocal Best Hits

```python
def reciprocal_best_hit(seq1, seq2, db="nr", program="blastp"):
    """Check if two sequences are reciprocal best BLAST hits."""
    from Bio.Blast import NCBIWWW, NCBIXML

    # Forward search
    handle1 = NCBIWWW.qblast(program, db, seq1)
    record1 = NCBIXML.read(handle1)
    handle1.close()

    if not record1.alignments:
        return False
    best_hit1 = record1.alignments[0].accession

    # Reverse search
    handle2 = NCBIWWW.qblast(program, db, seq2)
    record2 = NCBIXML.read(handle2)
    handle2.close()

    if not record2.alignments:
        return False
    best_hit2 = record2.alignments[0].accession

    return best_hit1 == seq2 and best_hit2 == seq1
```

## Error Handling

```python
from Bio.Blast import NCBIWWW, NCBIXML
from urllib.error import HTTPError
import time

def robust_blast(sequence, program="blastp", database="nr", max_retries=3):
    """Execute BLAST with retry logic."""
    for attempt in range(max_retries):
        try:
            handle = NCBIWWW.qblast(program, database, sequence)
            result = NCBIXML.read(handle)
            handle.close()
            return result
        except HTTPError as e:
            print(f"Attempt {attempt + 1} failed: HTTP {e.code}")
            if attempt < max_retries - 1:
                time.sleep(5 * (attempt + 1))
        except Exception as e:
            print(f"Error: {e}")
            raise

    raise RuntimeError("BLAST failed after maximum retries")
```

## Recommendations

1. **Use XML format** (outfmt 5) for reliable parsing
2. **Save results locally** - don't repeat searches unnecessarily
3. **Set reasonable E-value thresholds** (0.001-0.01 for most uses)
4. **Use local BLAST** for large-scale or repeated searches
5. **Filter by identity**, not just E-value
6. **Cache results** to avoid redundant API calls
7. **Handle empty results** gracefully
8. **Batch queries** when possible
9. **Consider DIAMOND** for very large protein searches
10. **Respect NCBI rate limits** for remote searches
