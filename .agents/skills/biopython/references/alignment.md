# Sequence Alignment

## Bio.Align: Pairwise Alignment

### PairwiseAligner Configuration

The `PairwiseAligner` class automatically selects the optimal algorithm (Needleman-Wunsch, Smith-Waterman, Gotoh, or Waterman-Smith-Beyer) based on scoring parameters.

```python
from Bio import Align

aligner = Align.PairwiseAligner()

# Default parameters (Biopython 1.85+):
# match: +1.0, mismatch: 0.0, gaps: -1.0
```

### Scoring Configuration

```python
# Simple scoring
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.gap_score = -2

# Affine gap penalties
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

# Different internal vs terminal gaps
aligner.internal_open_gap_score = -5
aligner.internal_extend_gap_score = -1

# Free end gaps (semi-global alignment)
aligner.left_open_gap_score = 0
aligner.left_extend_gap_score = 0
aligner.right_open_gap_score = 0
aligner.right_extend_gap_score = 0
```

### Alignment Modes

| Mode | Use Case |
|------|----------|
| `global` | Full-length comparison of similar sequences |
| `local` | Finding conserved regions within longer sequences |

```python
aligner.mode = 'local'  # Smith-Waterman
aligner.mode = 'global'  # Needleman-Wunsch (default)
```

### Executing Alignments

```python
from Bio.Seq import Seq

target = Seq("GAATTC")
query = Seq("GATTC")

# Get all optimal alignments
alignments = aligner.align(target, query)

for aln in alignments:
    print(aln)
    print(f"Score: {aln.score}")

# Score only (faster)
score = aligner.score(target, query)
```

### Substitution Matrices

For protein alignments, use established scoring matrices:

```python
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")
aligner.substitution_matrix = blosum62

# Available matrices
print(substitution_matrices.load())
```

| Matrix | Typical Use |
|--------|-------------|
| BLOSUM62 | General protein alignment |
| BLOSUM80 | Closely related sequences |
| BLOSUM45 | Distantly related sequences |
| PAM250 | Divergent proteins |
| PAM30 | Closely related proteins |

## Bio.AlignIO: Multiple Sequence Alignments

### Reading Alignments

```python
from Bio import AlignIO

# Single alignment
msa = AlignIO.read("cytochrome.aln", "clustal")

# Multiple alignments in one file
for alignment in AlignIO.parse("multi.aln", "clustal"):
    print(f"{len(alignment)} sequences, {alignment.get_alignment_length()} positions")
```

### Supported Formats

| Format | String | Notes |
|--------|--------|-------|
| Clustal | `"clustal"` | ClustalW/Omega output |
| PHYLIP | `"phylip"` | Strict 10-char names |
| PHYLIP (relaxed) | `"phylip-relaxed"` | Longer names allowed |
| Stockholm | `"stockholm"` | Pfam/Rfam format |
| FASTA | `"fasta"` | Aligned sequences with gaps |
| NEXUS | `"nexus"` | PAUP/MrBayes format |
| MSF | `"msf"` | GCG format |
| MAF | `"maf"` | Genomic alignments |

### Writing and Converting

```python
# Write alignment
AlignIO.write(msa, "output.phy", "phylip")

# Direct conversion
AlignIO.convert("input.aln", "clustal", "output.fasta", "fasta")
```

### Alignment Object Operations

```python
msa = AlignIO.read("proteins.aln", "clustal")

# Properties
print(f"Sequences: {len(msa)}")
print(f"Length: {msa.get_alignment_length()}")

# Access sequences
for record in msa:
    print(f"{record.id}: {record.seq[:50]}...")

# Column extraction
first_column = msa[:, 0]

# Region extraction
region = msa[:, 100:200]

# Single sequence
first_seq = msa[0]
```

### Alignment Statistics

```python
from Bio.Align import AlignInfo

summary = AlignInfo.SummaryInfo(msa)
consensus = summary.gap_consensus(threshold=0.7)
```

## Programmatic Alignment Creation

```python
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

aligned_records = [
    SeqRecord(Seq("MKTL--VELL"), id="seq1"),
    SeqRecord(Seq("MKTLIVELL-"), id="seq2"),
    SeqRecord(Seq("MKT--VELLL"), id="seq3"),
]

msa = MultipleSeqAlignment(aligned_records)
```

## External Aligners

### Clustal Omega

```python
from Bio.Align.Applications import ClustalOmegaCommandline

cmd = ClustalOmegaCommandline(
    infile="unaligned.fasta",
    outfile="aligned.aln",
    auto=True
)
stdout, stderr = cmd()
result = AlignIO.read("aligned.aln", "clustal")
```

### MUSCLE

```python
from Bio.Align.Applications import MuscleCommandline

cmd = MuscleCommandline(input="input.fasta", out="aligned.fasta")
stdout, stderr = cmd()
```

## Analysis Functions

### Pairwise Identity Calculation

```python
def sequence_identity(seq1, seq2):
    """Calculate percent identity between aligned sequences."""
    paired = [(a, b) for a, b in zip(seq1, seq2) if a != '-' and b != '-']
    if not paired:
        return 0.0
    matches = sum(1 for a, b in paired if a == b)
    return matches / len(paired)

# Calculate all pairwise identities
for i, rec1 in enumerate(msa):
    for rec2 in msa[i+1:]:
        pct = sequence_identity(rec1.seq, rec2.seq)
        print(f"{rec1.id} vs {rec2.id}: {pct:.1%}")
```

### Identify Conserved Positions

```python
def find_conserved_columns(msa, threshold=0.9):
    """Find columns with conservation above threshold."""
    conserved = []
    for col_idx in range(msa.get_alignment_length()):
        column = msa[:, col_idx]
        residue_counts = {}
        for residue in column:
            residue_counts[residue] = residue_counts.get(residue, 0) + 1
        max_freq = max(residue_counts.values()) / len(column)
        if max_freq >= threshold:
            conserved.append(col_idx)
    return conserved

conserved_cols = find_conserved_columns(msa, 0.8)
```

### Remove Gap-Only Columns

```python
def remove_empty_columns(msa):
    """Remove columns that are entirely gaps."""
    valid_cols = []
    for i in range(msa.get_alignment_length()):
        col = msa[:, i]
        if set(col) != {'-'}:
            valid_cols.append(i)

    # Extract valid columns
    new_seqs = []
    for record in msa:
        new_seq = ''.join(str(record.seq)[i] for i in valid_cols)
        new_seqs.append(SeqRecord(Seq(new_seq), id=record.id))

    return MultipleSeqAlignment(new_seqs)
```

## Examples

### Local Alignment for Domain Detection

```python
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

aligner = PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 3
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

full_protein = Seq("MKTLLVELLKRGGGGMKTLLVELLKRGGGG")
domain = Seq("MKTLLVELLKR")

alignments = aligner.align(full_protein, domain)
best = alignments[0]
print(best)
```

### Protein Alignment with BLOSUM62

```python
from Bio.Align import PairwiseAligner, substitution_matrices

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

human = Seq("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
mouse = Seq("MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH")

alignments = aligner.align(human, mouse)
print(alignments[0])
print(f"Score: {alignments[0].score}")
```

## Recommendations

1. Select appropriate scoring parameters for your sequences
2. Use `global` for comparing full-length homologs, `local` for finding conserved domains
3. Apply BLOSUM62 for general protein work; BLOSUM80 for close homologs
4. Higher gap penalties produce fewer, longer gaps
5. Validate alignment quality manually for important analyses
6. Use external tools (ClustalO, MUSCLE) for large MSAs
