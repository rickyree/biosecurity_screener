# scikit-bio API Reference

Detailed method signatures, extended examples, and troubleshooting.

## Contents

1. [Sequence Classes](#sequence-classes)
2. [Alignment](#alignment)
3. [Phylogenetic Trees](#phylogenetic-trees)
4. [Diversity Metrics](#diversity-metrics)
5. [Ordination](#ordination)
6. [Statistical Tests](#statistical-tests)
7. [Distance Matrices](#distance-matrices)
8. [File I/O](#file-io)
9. [Common Issues](#common-issues)

---

## Sequence Classes

### Creating Sequences

```python
from skbio import DNA, RNA, Protein, Sequence

# With metadata
gene = DNA('ATGCGATCGATCG', metadata={'id': 'gene_001', 'organism': 'E. coli'})
transcript = RNA('AUGCGAUCGAUCG')
enzyme = Protein('MKTAYIAKQRQISFVK')

# Transformations
antisense = gene.reverse_complement()
mrna = gene.transcribe()
peptide = mrna.translate()

# Alternate genetic codes
bacterial_peptide = mrna.translate(genetic_code=11)
```

### Pattern Matching

```python
import re

seq = DNA('ATGCATGCATGCATGC')

# Regex-based search
start_codons = seq.find_with_regex('ATG.{6}')

# Manual iteration
for match in re.finditer('ATG', str(seq)):
    print(f"Position: {match.start()}")

# K-mer frequencies
kmer_counts = seq.kmer_frequencies(k=4)
```

### Metadata Handling

```python
# Sequence-level
seq = DNA('ATCG', metadata={'id': 'contig_17', 'length': 4})
print(seq.metadata['id'])

# Positional (quality scores from FASTQ)
fastq_seq = DNA.read('illumina.fastq', format='fastq', phred_offset=33)
scores = fastq_seq.positional_metadata['quality']

# Interval (annotations)
seq.interval_metadata.add([(10, 50)], metadata={'feature': 'promoter'})
```

### Distance Between Sequences

```python
from skbio import DNA
from skbio.sequence.distance import kmer_distance

seqA = DNA('ATCGATCGATCG')
seqB = DNA('ATCG--CGATCG')

# Hamming distance
hamming = seqA.distance(seqB)

# K-mer based
kmer_dist = seqA.distance(seqB, metric=kmer_distance)
```

---

## Alignment

### Pairwise Alignment

```python
from skbio.alignment import local_pairwise_align_ssw, global_pairwise_align
from skbio import DNA, Protein

# Local (SSW-accelerated)
query = DNA('ATCGATCGATCGATCG')
target = DNA('ATCGAAAATCGATCG')
local_aln = local_pairwise_align_ssw(query, target)

# Inspect results
print(f"Score: {local_aln.score}")
print(f"Query start: {local_aln.target_begin}")

# Global with custom scoring
from skbio.alignment import AlignScorer

scorer = AlignScorer(
    match_score=3,
    mismatch_score=-2,
    gap_open_penalty=6,
    gap_extend_penalty=1
)
global_aln = global_pairwise_align(query, target, scorer=scorer)
```

### Protein Alignment with Substitution Matrix

```python
from skbio.alignment import StripedSmithWaterman
from skbio import Protein

prot_query = Protein('MVLSPADKTNVKAAW')
prot_target = Protein('MVLSGEDKSNIKAAW')

aligner = StripedSmithWaterman(
    str(prot_query),
    gap_open_penalty=11,
    gap_extend_penalty=1,
    substitution_matrix='blosum62'
)
prot_aln = aligner(str(prot_target))
```

### Multiple Sequence Alignment

```python
from skbio.alignment import TabularMSA
from skbio import DNA

# Read from file
msa = TabularMSA.read('aligned_genes.fasta', constructor=DNA)

# Build manually
aligned = TabularMSA([
    DNA('ATCG--AA'),
    DNA('ATGG--AA'),
    DNA('ATCGATAA')
])

# Operations
consensus = msa.consensus()
majority = msa.majority_consensus()

# Access data
first_seq = msa[0]
third_column = msa[:, 2]

# Remove gappy columns
filtered = msa.omit_gap_positions(maximum_gap_frequency=0.4)

# Per-position entropy
entropies = msa.position_entropies()
```

### CIGAR Strings

```python
from skbio.alignment import AlignPath

# Parse CIGAR
cigar_str = "12M3I8M2D15M"
path = AlignPath.from_cigar(cigar_str, target_length=120, query_length=80)

# Generate CIGAR from alignment
output_cigar = local_aln.to_cigar()
```

---

## Phylogenetic Trees

### Construction

```python
from skbio import DistanceMatrix
from skbio.tree import nj, upgma, bme

dm = DistanceMatrix([
    [0, 3, 8, 8],
    [3, 0, 9, 9],
    [8, 9, 0, 6],
    [8, 9, 6, 0]
], ids=['alpha', 'beta', 'gamma', 'delta'])

# Neighbor joining
tree_nj = nj(dm)

# UPGMA (molecular clock)
tree_upgma = upgma(dm)

# Balanced minimum evolution (large datasets)
tree_bme = bme(dm)
```

### Manipulation

```python
from skbio import TreeNode

tree = TreeNode.read('mammals.nwk', format='newick')

# Traversal modes
for node in tree.preorder():
    if node.name:
        print(node.name)

# Tips only
tip_names = [t.name for t in tree.tips()]

# Locate node
target = tree.find('homo_sapiens')

# Midpoint rooting
rooted = tree.root_at_midpoint()

# Subset to specific taxa
reduced = tree.shear(['mus_musculus', 'rattus_norvegicus', 'homo_sapiens'])

# Most recent common ancestor
mrca = tree.lowest_common_ancestor(['mus_musculus', 'rattus_norvegicus'])
subtree = mrca.copy()
```

### Modifying Trees

```python
# Add node
parent = tree.find('rodentia')
new_child = TreeNode(name='new_species', length=0.05)
parent.append(new_child)

# Remove node
obsolete = tree.find('deprecated_taxon')
obsolete.parent.remove(obsolete)
```

### Distance Calculations

```python
# Path length between nodes
node_a = tree.find('species_a')
node_b = tree.find('species_b')
path_dist = node_a.distance(node_b)

# All pairwise distances
all_pairs = tree.cophenetic_matrix()

# Topology comparison
rf_dist = tree.robinson_foulds(reference_tree)
```

### Visualization

```python
# Text representation
print(tree.ascii_art())

# Export for external tools
tree.write('output.nwk', format='newick')
```

---

## Diversity Metrics

### Alpha Diversity

```python
from skbio.diversity import alpha_diversity, get_alpha_diversity_metrics
import numpy as np

# Abundance data
counts = np.array([
    [25, 10, 0, 5],
    [8, 0, 15, 12],
    [12, 12, 8, 8]
])
sample_names = ['sampleA', 'sampleB', 'sampleC']

# Available metrics
print(get_alpha_diversity_metrics())

# Calculate
shannon = alpha_diversity('shannon', counts, ids=sample_names)
simpson = alpha_diversity('simpson', counts, ids=sample_names)
richness = alpha_diversity('observed_otus', counts, ids=sample_names)
chao1_est = alpha_diversity('chao1', counts, ids=sample_names)
```

### Phylogenetic Alpha Diversity

```python
from skbio import TreeNode

phylo = TreeNode.read('otus.nwk')
otu_names = ['otu_1', 'otu_2', 'otu_3', 'otu_4']

faith = alpha_diversity('faith_pd', counts, ids=sample_names,
                        tree=phylo, otu_ids=otu_names)
```

### Beta Diversity

```python
from skbio.diversity import beta_diversity, partial_beta_diversity

# Standard metrics
bray = beta_diversity('braycurtis', counts, ids=sample_names)
jaccard = beta_diversity('jaccard', counts, ids=sample_names)

# Phylogenetic metrics
unweighted_uf = beta_diversity('unweighted_unifrac', counts,
                                ids=sample_names,
                                tree=phylo,
                                otu_ids=otu_names)

weighted_uf = beta_diversity('weighted_unifrac', counts,
                              ids=sample_names,
                              tree=phylo,
                              otu_ids=otu_names)
```

### Selective Computation

```python
# Only specific pairs
pairs_of_interest = [('sampleA', 'sampleB'), ('sampleA', 'sampleC')]
partial_dm = partial_beta_diversity('braycurtis', counts,
                                    ids=sample_names,
                                    id_pairs=pairs_of_interest)
```

### Rarefaction

```python
from skbio.diversity import subsample_counts
import numpy as np

# Rarefy each sample to fixed depth
target_depth = 500
rarefied = np.array([subsample_counts(row, n=target_depth) for row in counts])

# Repeated rarefaction for confidence intervals
iterations = 100
shannon_replicates = []
for _ in range(iterations):
    rare_counts = np.array([subsample_counts(row, n=target_depth) for row in counts])
    shannon_replicates.append(alpha_diversity('shannon', rare_counts))

mean_shannon = np.mean(shannon_replicates, axis=0)
sd_shannon = np.std(shannon_replicates, axis=0)
```

---

## Ordination

### PCoA

```python
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix

dm = DistanceMatrix(...)
result = pcoa(dm)

# Coordinates
x = result.samples['PC1']
y = result.samples['PC2']

# Variance explained
pct_var = result.proportion_explained

# Eigenvalues
eig = result.eigvals

# Persistence
result.write('pcoa_output.txt')
```

### Plotting PCoA

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(x, y)
ax.set_xlabel(f'PC1 ({pct_var[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pct_var[1]*100:.1f}%)')
plt.tight_layout()
```

### CCA

```python
from skbio.stats.ordination import cca
import pandas as pd
import numpy as np

# Species abundances
species = np.array([
    [15, 8, 2],
    [3, 12, 7],
    [9, 9, 6]
])

# Environmental gradients
env = pd.DataFrame({
    'moisture': [35, 60, 48],
    'temp_c': [18, 24, 21],
    'elevation': [200, 350, 275]
})

result = cca(species, env,
             sample_ids=['loc_1', 'loc_2', 'loc_3'],
             species_ids=['sp_X', 'sp_Y', 'sp_Z'])

cca1 = result.samples['CCA1']
cca2 = result.samples['CCA2']
env_arrows = result.biplot_scores
```

### RDA

```python
from skbio.stats.ordination import rda

# Linear model variant
rda_result = rda(species, env,
                 sample_ids=['loc_1', 'loc_2', 'loc_3'],
                 species_ids=['sp_X', 'sp_Y', 'sp_Z'])
```

---

## Statistical Tests

### PERMANOVA

```python
from skbio.stats.distance import permanova

dm = DistanceMatrix(...)
groups = ['treatment', 'treatment', 'control', 'control', 'treatment', 'control']

result = permanova(dm, groups, permutations=999)

print(f"F-statistic: {result['test statistic']:.3f}")
print(f"p-value: {result['p-value']:.4f}")
print(f"n: {result['sample size']}")
print(f"Groups: {result['number of groups']}")
```

### ANOSIM

```python
from skbio.stats.distance import anosim

result = anosim(dm, groups, permutations=999)
print(f"R: {result['test statistic']:.3f}")
print(f"p: {result['p-value']:.4f}")
```

### PERMDISP

```python
from skbio.stats.distance import permdisp

# Homogeneity of group dispersions
result = permdisp(dm, groups, permutations=999)
print(f"F: {result['test statistic']:.3f}")
print(f"p: {result['p-value']:.4f}")
```

### Mantel Test

```python
from skbio.stats.distance import mantel

dm_genetic = DistanceMatrix(...)
dm_geographic = DistanceMatrix(...)

r, p, n = mantel(dm_genetic, dm_geographic, method='pearson', permutations=999)
print(f"Correlation: {r:.3f}")
print(f"p-value: {p:.4f}")

# Spearman variant
r_sp, p_sp, n_sp = mantel(dm_genetic, dm_geographic, method='spearman', permutations=999)
```

---

## Distance Matrices

### Construction

```python
from skbio import DistanceMatrix, DissimilarityMatrix
import numpy as np

# Symmetric distances
data = np.array([
    [0.0, 0.2, 0.5],
    [0.2, 0.0, 0.4],
    [0.5, 0.4, 0.0]
])
dm = DistanceMatrix(data, ids=['x', 'y', 'z'])

# Asymmetric dissimilarities
asym = np.array([
    [0.0, 0.3, 0.6],
    [0.4, 0.0, 0.5],
    [0.7, 0.6, 0.0]
])
dissim = DissimilarityMatrix(asym, ids=['p', 'q', 'r'])
```

### Access and Manipulation

```python
# Element access
val = dm['x', 'y']
row = dm['x']

# Subset
subset = dm.filter(['x', 'z'])

# I/O
dm.write('pairwise.txt')
loaded = DistanceMatrix.read('pairwise.txt')

# Scipy format
condensed = dm.condensed_form()

# Pandas
df = dm.to_data_frame()
```

---

## File I/O

### Reading Sequences

```python
import skbio

# Single sequence
seq = skbio.DNA.read('reference.fasta', format='fasta')

# Iterator for large files
for record in skbio.io.read('reads.fasta', format='fasta', constructor=skbio.DNA):
    print(record.metadata['id'], len(record))

# Collect all
all_seqs = list(skbio.io.read('genes.fasta', format='fasta', constructor=skbio.DNA))
```

### Reading FASTQ with Quality

```python
for seq in skbio.io.read('hiseq.fastq', format='fastq', constructor=skbio.DNA):
    qual = seq.positional_metadata['quality']
    avg_q = qual.mean()
    print(f"{seq.metadata['id']}: mean Q={avg_q:.1f}")
```

### Writing Sequences

```python
# Single
seq.write('gene.fasta', format='fasta')

# Multiple
skbio.io.write([seq1, seq2, seq3], format='fasta', into='combined.fasta')

# Custom line width
seq.write('formatted.fasta', format='fasta', max_width=80)
```

### BIOM Tables

```python
from skbio import Table

# Load
tbl = Table.read('features.biom', format='hdf5')

# Inspect
samples = tbl.ids(axis='sample')
features = tbl.ids(axis='observation')
matrix = tbl.matrix_data.toarray()

# Filter samples
abundant = tbl.filter(lambda row, id_, md: row.sum() > 1000, axis='sample')

# Filter features
prevalent = tbl.filter(lambda col, id_, md: (col > 0).sum() >= 5,
                       axis='observation')

# Normalize
rel_abund = tbl.norm(axis='sample', inplace=False)

# Save
tbl.write('processed.biom', format='hdf5')
```

### Format Conversion

```python
# FASTQ to FASTA
seqs = skbio.io.read('raw.fastq', format='fastq', constructor=skbio.DNA)
skbio.io.write(seqs, format='fasta', into='converted.fasta')

# GenBank to FASTA
gb_seqs = skbio.io.read('annotations.gb', format='genbank', constructor=skbio.DNA)
skbio.io.write(gb_seqs, format='fasta', into='sequences.fasta')
```

---

## Common Issues

### Duplicate IDs

```python
# Error: ValueError: Ids must be unique
# Fix: Deduplicate

seen = set()
unique = []
for seq in sequences:
    seq_id = seq.metadata['id']
    if seq_id not in seen:
        unique.append(seq)
        seen.add(seq_id)
```

### Non-integer Counts

```python
# Error: ValueError: Counts must be integers
# Fix: Convert proportions to counts

int_counts = (proportions * 10000).astype(int)
```

### Memory with Large Files

```python
# Problem: MemoryError
# Fix: Process iteratively

for seq in skbio.io.read('huge.fasta', format='fasta', constructor=skbio.DNA):
    result = process_sequence(seq)
    save_result(result)
```

### Tree-Table ID Mismatch

```python
# Problem: Tips don't match feature IDs
# Fix: Verify and align

tree_tips = {t.name for t in tree.tips()}
table_features = set(feature_ids)

missing_in_tree = table_features - tree_tips
missing_in_table = tree_tips - table_features

# Prune tree to match table
tree_matched = tree.shear(list(table_features))
```

### Alignment Length Mismatch

```python
# Problem: Sequences have different lengths (pre-aligned)
# Fix: Remove gaps first

clean_a = seq_a.degap()
clean_b = seq_b.degap()
alignment = local_pairwise_align_ssw(clean_a, clean_b)
```

---

## Performance Guidelines

| Scenario | Recommendation |
|----------|----------------|
| Large sequence files | Use generators instead of `list()` |
| Big BIOM tables | Prefer HDF5 over JSON format |
| Partial comparisons | Use `partial_beta_diversity()` |
| Huge phylogenies | Use BME instead of NJ |
| Exploratory analysis | Subsample data first |
| Repeated analyses | Cache distance matrices and ordinations |

---

## Integration Patterns

### With pandas

```python
import pandas as pd
from skbio import DistanceMatrix

dm = DistanceMatrix(...)
df = dm.to_data_frame()

alpha = alpha_diversity('shannon', counts, ids=sample_ids)
alpha_df = pd.DataFrame({'shannon_index': alpha})
```

### With matplotlib/seaborn

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Ordination scatter
fig, ax = plt.subplots()
scatter = ax.scatter(pc1, pc2, c=group_labels, cmap='Set2')
ax.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}%)')
plt.colorbar(scatter)

# Distance heatmap
sns.heatmap(dm.to_data_frame(), cmap='YlOrRd', square=True)
```

### With QIIME 2

```python
# Export from QIIME 2
# qiime tools export --input-path data.qza --output-path exported/

# Process with scikit-bio
tbl = Table.read('exported/feature-table.biom')
# ... analysis ...
tbl.write('results.biom')

# Import back
# qiime tools import --input-path results.biom --output-path results.qza
```
