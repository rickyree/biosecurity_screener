# Advanced Biopython Features

## Sequence Motifs (Bio.motifs)

### Creating Motifs from Sequences

```python
from Bio import motifs
from Bio.Seq import Seq

instances = [
    Seq("TACGC"),
    Seq("TACCC"),
    Seq("TATGC"),
    Seq("AACGC"),
    Seq("TACGA"),
]

motif = motifs.create(instances)
```

### Consensus Sequences

```python
# Standard consensus
print(motif.counts.consensus)

# Degenerate consensus (IUPAC ambiguity)
print(motif.counts.degenerate_consensus)
```

### Position Weight Matrix

```python
# Normalize with pseudocounts
pwm = motif.counts.normalize(pseudocounts=0.5)

# Position-specific scoring matrix
pssm = pwm.log_odds()

# Information content
ic = motif.counts.information_content()
```

### Scanning Sequences

```python
target = Seq("ATCGTACGCATTATGCCCATACGAATTT")

for position, score in pssm.search(target, threshold=3.0):
    print(f"Position {position}: score {score:.2f}")
```

### Reading Motif Files

```python
# JASPAR format
with open("motif.jaspar") as f:
    m = motifs.read(f, "jaspar")

# Parse multiple motifs
with open("database.jaspar") as f:
    for m in motifs.parse(f, "jaspar"):
        print(m.name)
```

## Sequence Utilities (Bio.SeqUtils)

### GC Content

```python
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq

seq = Seq("ATCGATCGATCGATCG")
gc = gc_fraction(seq)
print(f"GC content: {gc:.1%}")
```

### Molecular Weight

```python
from Bio.SeqUtils import molecular_weight

# DNA
dna_mw = molecular_weight(Seq("ATCGATCG"), seq_type="DNA")

# Protein
protein_mw = molecular_weight(Seq("MKFLILFFVLFFVLFFAVQALQDL"), seq_type="protein")
```

### Melting Temperature

```python
from Bio.SeqUtils import MeltingTemp as mt

primer = Seq("GCTAGCTAGCTAGCTA")

# Nearest-neighbor method
tm_nn = mt.Tm_NN(primer)

# With salt correction
tm_salt = mt.Tm_NN(primer, Na=50, Mg=1.5)

# Wallace rule (rough estimate)
tm_wallace = mt.Tm_Wallace(primer)
```

### Protein Analysis

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis

protein = "MKFLILFFVLFFVLFFAVQALQDLMRIIKEN"
analysis = ProteinAnalysis(protein)

print(f"MW: {analysis.molecular_weight():.1f} Da")
print(f"pI: {analysis.isoelectric_point():.2f}")
print(f"Instability: {analysis.instability_index():.1f}")
print(f"Aromaticity: {analysis.aromaticity():.3f}")
print(f"GRAVY: {analysis.gravy():.3f}")

# Secondary structure propensity
helix, turn, sheet = analysis.secondary_structure_fraction()
print(f"Helix: {helix:.1%}, Turn: {turn:.1%}, Sheet: {sheet:.1%}")

# Amino acid composition
composition = analysis.get_amino_acids_percent()

# Extinction coefficient (reduced cysteines)
ext_coef = analysis.molar_extinction_coefficient()
```

## Restriction Enzymes (Bio.Restriction)

### Finding Cut Sites

```python
from Bio import Restriction
from Bio.Seq import Seq

dna = Seq("GAATTCATCGATCGATGAATTCGGATCCAAGCTT")

# Single enzyme
ecori_sites = Restriction.EcoRI.search(dna)
print(f"EcoRI cuts at: {ecori_sites}")

# Multiple enzymes
batch = Restriction.RestrictionBatch(["EcoRI", "BamHI", "HindIII"])
results = batch.search(dna)

for enzyme, sites in results.items():
    if sites:
        print(f"{enzyme}: {sites}")
```

### Finding Enzymes That Cut

```python
# Which enzymes cut this sequence?
rb = Restriction.CommOnly  # Commonly available enzymes
analysis = Restriction.Analysis(rb, dna)

# Enzymes that cut
cutters = analysis.with_sites()

# Enzymes that don't cut
non_cutters = analysis.without_sites()
```

## Genetic Codes (Bio.Data.CodonTable)

```python
from Bio.Data import CodonTable

# Standard genetic code (table 1)
standard = CodonTable.unambiguous_dna_by_id[1]
print(standard)

# Mitochondrial code (table 2)
mito = CodonTable.unambiguous_dna_by_id[2]

# Look up codon
amino_acid = standard.forward_table["ATG"]

# Start and stop codons
starts = standard.start_codons
stops = standard.stop_codons
```

## Population Genetics (Bio.PopGen)

### GenePop Files

```python
from Bio.PopGen import GenePop

with open("population.gen") as f:
    record = GenePop.read(f)

print(f"Populations: {len(record.populations)}")
print(f"Loci: {record.loci_list}")

for i, pop in enumerate(record.populations):
    print(f"Population {i+1}: {len(pop)} individuals")
```

### Population Statistics

```python
from Bio.PopGen.GenePop.Controller import GenePopController

ctrl = GenePopController()

# Allele frequencies
freqs = ctrl.calc_allele_genotype_freqs("data.gen")

# Fst calculation
fst = ctrl.calc_fst_all("data.gen")

# Hardy-Weinberg test
hw = ctrl.test_hw_pop("data.gen", "probability")
```

## Sequence Features (Bio.SeqFeature)

### Creating Features

```python
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Simple feature
cds = SeqFeature(
    location=FeatureLocation(start=100, end=400),
    type="CDS",
    strand=1,
    qualifiers={"gene": ["myGene"], "product": ["hypothetical protein"]}
)

# Add to record
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

record = SeqRecord(Seq("N" * 500), id="seq001")
record.features.append(cds)
```

### Extracting Feature Sequences

```python
# Extract sequence from feature location
feature_seq = cds.extract(record.seq)
```

### Complex Locations

```python
from Bio.SeqFeature import CompoundLocation

# Join (spliced feature)
exon1 = FeatureLocation(100, 200)
exon2 = FeatureLocation(300, 400)
spliced = CompoundLocation([exon1, exon2])

feature = SeqFeature(location=spliced, type="mRNA")
```

## Ambiguity Codes

```python
from Bio.Data import IUPACData

# DNA ambiguity codes
print(IUPACData.ambiguous_dna_letters)  # ACGTRYSWKMBDHVN

# What does R mean?
print(IUPACData.ambiguous_dna_values["R"])  # ['A', 'G']

# Protein ambiguity
print(IUPACData.ambiguous_protein_letters)
```

## FASTQ Quality Scores

```python
from Bio import SeqIO

for record in SeqIO.parse("reads.fastq", "fastq"):
    quals = record.letter_annotations["phred_quality"]

    avg_qual = sum(quals) / len(quals)
    min_qual = min(quals)

    print(f"{record.id}: avg={avg_qual:.1f}, min={min_qual}")

    # Quality filtering
    if min_qual >= 20:
        # Process high-quality read
        pass
```

## Genome Diagrams (GenomeDiagram)

```python
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from reportlab.lib import colors

# Load annotated sequence
record = SeqIO.read("genome.gb", "genbank")

# Create diagram
diagram = GenomeDiagram.Diagram("Genome Map")
track = diagram.new_track(1, greytrack=True)
feature_set = track.new_set()

# Add features with colors
for feature in record.features:
    if feature.type == "CDS":
        color = colors.blue
    elif feature.type == "gene":
        color = colors.lightblue
    else:
        continue

    feature_set.add_feature(
        feature,
        color=color,
        label=True,
        label_size=8
    )

# Render
diagram.draw(format="linear", pagesize="A4", fragments=1)
diagram.write("genome_map.pdf", "PDF")
```

## Clustering (Bio.Cluster)

```python
from Bio.Cluster import kcluster
import numpy as np

# Expression data matrix (genes x conditions)
data = np.array([
    [1.2, 0.8, 2.1],
    [1.1, 0.9, 2.0],
    [0.2, 2.5, 0.3],
    [0.3, 2.4, 0.2],
])

# K-means clustering
cluster_ids, error, n_found = kcluster(data, nclusters=2)
print(f"Cluster assignments: {cluster_ids}")
```

## Practical Examples

### Find Open Reading Frames

```python
from Bio.Seq import Seq

def find_orfs(sequence, min_protein_length=50):
    """Find all ORFs in a sequence."""
    orfs = []

    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate()

            start = 0
            while start < len(trans):
                # Find start codon position
                met_pos = str(trans).find("M", start)
                if met_pos == -1:
                    break

                # Find stop codon
                stop_pos = str(trans).find("*", met_pos)
                if stop_pos == -1:
                    stop_pos = len(trans)

                protein_len = stop_pos - met_pos
                if protein_len >= min_protein_length:
                    nt_start = frame + met_pos * 3
                    nt_end = frame + stop_pos * 3

                    orfs.append({
                        'strand': strand,
                        'frame': frame,
                        'start': nt_start,
                        'end': nt_end,
                        'length': protein_len
                    })

                start = stop_pos + 1

    return orfs
```

### Codon Usage Analysis

```python
from Bio import SeqIO
from collections import Counter

def codon_usage(fasta_file):
    """Calculate codon frequencies from coding sequences."""
    codons = Counter()

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        # Trim to multiple of 3
        seq = seq[:len(seq) - len(seq) % 3]

        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if "N" not in codon:  # Skip ambiguous
                codons[codon] += 1

    total = sum(codons.values())
    return {codon: count / total for codon, count in codons.items()}
```

### Sequence Complexity (Shannon Entropy)

```python
import math
from collections import Counter

def kmer_complexity(sequence, k=2):
    """Calculate sequence complexity via k-mer entropy."""
    seq = str(sequence)
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]

    counts = Counter(kmers)
    total = len(kmers)

    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)

    # Normalize by max possible entropy
    max_entropy = math.log2(4 ** k)

    return entropy / max_entropy if max_entropy > 0 else 0
```

### Promoter Region Extraction

```python
from Bio import SeqIO

def extract_promoters(genbank_file, upstream=500):
    """Extract upstream regions of annotated genes."""
    record = SeqIO.read(genbank_file, "genbank")
    promoters = []

    for feature in record.features:
        if feature.type != "gene":
            continue

        gene_name = feature.qualifiers.get("gene", ["unknown"])[0]

        if feature.strand == 1:
            start = max(0, feature.location.start - upstream)
            end = feature.location.start
            prom_seq = record.seq[start:end]
        else:
            start = feature.location.end
            end = min(len(record), feature.location.end + upstream)
            prom_seq = record.seq[start:end].reverse_complement()

        promoters.append({
            'gene': gene_name,
            'sequence': prom_seq
        })

    return promoters
```

## Recommendations

1. **Combine modules** for complex workflows
2. **Use pseudocounts** in motif analysis to avoid zero probabilities
3. **Validate inputs** - check file formats and data quality
4. **Handle edge cases** - empty results, missing data
5. **Cache expensive computations** for large analyses
6. **Use appropriate genetic codes** for non-standard organisms
7. **Document parameters** for reproducibility
8. **Test with known data** before production runs
