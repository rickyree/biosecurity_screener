# Sequence Objects and File I/O

## Bio.Seq: The Sequence Object

### Creating Sequences

```python
from Bio.Seq import Seq

nucleotide = Seq("GATTACA")
print(len(nucleotide))      # 7
print(nucleotide[2:5])      # TTA
print(nucleotide.upper())   # GATTACA
```

### Biological Operations

| Method | Description | Example |
|--------|-------------|---------|
| `complement()` | Complementary strand | `GATTACA` → `CTAATGT` |
| `reverse_complement()` | Reverse complement | `GATTACA` → `TGTAATC` |
| `transcribe()` | DNA → RNA | `GATTACA` → `GAUUACA` |
| `back_transcribe()` | RNA → DNA | `GAUUACA` → `GATTACA` |
| `translate()` | Nucleotide → Protein | `ATGAAA` → `MK` |
| `translate(table=N)` | Use genetic code N | Mitochondrial: `table=2` |
| `translate(to_stop=True)` | Stop at first * | Truncate at stop codon |

```python
dna = Seq("ATGAAATTTGGATAG")
rna = dna.transcribe()
protein = dna.translate()
protein_truncated = dna.translate(to_stop=True)
```

## Bio.SeqIO: File Input/Output

### Core Functions

**`SeqIO.parse()`** - Iterate through multi-record files:

```python
from Bio import SeqIO

for record in SeqIO.parse("genome.fasta", "fasta"):
    print(f"{record.id}: {len(record.seq)} bp")
```

**`SeqIO.read()`** - Single-record files (validates exactly one record):

```python
single_record = SeqIO.read("plasmid.gb", "genbank")
```

**`SeqIO.write()`** - Write records to file:

```python
records = list(SeqIO.parse("input.fasta", "fasta"))
count = SeqIO.write(records, "output.fasta", "fasta")
```

**`SeqIO.convert()`** - Direct format conversion:

```python
SeqIO.convert("data.gbk", "genbank", "data.fasta", "fasta")
```

### Supported Formats

| Format | String | Notes |
|--------|--------|-------|
| FASTA | `"fasta"` | Most common sequence format |
| GenBank | `"genbank"` or `"gb"` | Rich annotation format |
| FASTQ | `"fastq"` | Includes quality scores |
| EMBL | `"embl"` | European database format |
| Swiss-Prot | `"swiss"` | Protein database format |
| Two-line FASTA | `"fasta-2line"` | Sequence on single line |
| Tab-delimited | `"tab"` | Simple ID-sequence pairs |

### SeqRecord Attributes

```python
record = SeqIO.read("gene.gb", "genbank")

record.id              # Primary identifier
record.name            # Short name
record.description     # Full description line
record.seq             # Seq object
record.annotations     # Dict: organism, date, etc.
record.features        # List of SeqFeature objects
record.letter_annotations  # Per-position data (quality scores)
```

### Modifying Records

```python
# Update attributes
record.id = "modified_id"
record.description = "Updated description"

# Slice preserves annotations
subset = record[100:500]

# Transform sequence
record.seq = record.seq.reverse_complement()
```

## Handling Large Files

### Memory-Efficient Iteration

```python
# Process one record at a time
for record in SeqIO.parse("large_genome.fasta", "fasta"):
    if some_condition(record):
        process(record)
```

### Random Access Strategies

**In-memory dictionary** (small files):

```python
record_dict = SeqIO.to_dict(SeqIO.parse("sequences.fasta", "fasta"))
target = record_dict["gene_001"]
```

**Indexed access** (medium files):

```python
idx = SeqIO.index("sequences.fasta", "fasta")
target = idx["gene_001"]
idx.close()
```

**SQLite-backed index** (very large files / multi-file):

```python
idx = SeqIO.index_db("cache.sqlite", "sequences.fasta", "fasta")
target = idx["gene_001"]
idx.close()
```

### High-Performance Parsers

For maximum speed with simple data extraction:

```python
from Bio.SeqIO.FastaIO import SimpleFastaParser

with open("reads.fasta") as handle:
    for title, sequence in SimpleFastaParser(handle):
        if len(sequence) > 1000:
            print(title)
```

```python
from Bio.SeqIO.QualityIO import FastqGeneralIterator

with open("reads.fastq") as handle:
    for title, seq, qual in FastqGeneralIterator(handle):
        avg_qual = sum(ord(c) - 33 for c in qual) / len(qual)
```

## Compressed Files

Automatic handling of gzip:

```python
for record in SeqIO.parse("genome.fasta.gz", "fasta"):
    print(record.id)
```

BGZF for random access to compressed files:

```python
from Bio import bgzf

with bgzf.open("data.fasta.bgz", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        process(record)
```

## Common Operations

### Filter by Criteria

```python
long_seqs = (r for r in SeqIO.parse("all.fasta", "fasta") if len(r) > 500)
SeqIO.write(long_seqs, "long_only.fasta", "fasta")
```

### Quality Filtering (FASTQ)

```python
high_quality = (
    rec for rec in SeqIO.parse("raw.fastq", "fastq")
    if min(rec.letter_annotations["phred_quality"]) >= 20
)
SeqIO.write(high_quality, "filtered.fastq", "fastq")
```

### Extract Metadata from GenBank

```python
for record in SeqIO.parse("collection.gbk", "genbank"):
    organism = record.annotations.get("organism", "Unknown")
    taxonomy = record.annotations.get("taxonomy", [])
    print(f"{record.id}: {organism}")
```

### Create Records Programmatically

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

new_record = SeqRecord(
    Seq("ATGCGATCGATCGATCG"),
    id="synthetic_001",
    name="SynGene",
    description="Synthetic test sequence"
)

SeqIO.write([new_record], "synthetic.fasta", "fasta")
```

### Batch Format Conversion

```python
input_files = ["sample1.gb", "sample2.gb", "sample3.gb"]
for infile in input_files:
    outfile = infile.replace(".gb", ".fasta")
    SeqIO.convert(infile, "genbank", outfile, "fasta")
```

## Recommendations

1. Use iterators (`SeqIO.parse`) rather than loading entire files
2. Choose `SeqIO.index()` for repeated random access
3. Use `SeqIO.index_db()` for millions of sequences
4. Apply low-level parsers for HTS data when speed matters
5. Cache fetched data locally instead of repeated downloads
6. Close indexed files explicitly or use context managers
7. Format strings are always lowercase: `"fasta"` not `"FASTA"`
