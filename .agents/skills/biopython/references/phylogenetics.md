# Phylogenetic Trees with Bio.Phylo

## Supported Formats

| Format | String | Description |
|--------|--------|-------------|
| Newick | `"newick"` | Standard parenthetical notation |
| NEXUS | `"nexus"` | Extended format with metadata |
| phyloXML | `"phyloxml"` | XML-based with rich annotations |
| NeXML | `"nexml"` | Modern XML format |
| CDAO | `"cdao"` | Ontology-based format |

## Reading and Writing Trees

### Load Trees

```python
from Bio import Phylo

# Single tree
tree = Phylo.read("species.nwk", "newick")

# Multiple trees
trees = list(Phylo.parse("bootstrap_trees.nwk", "newick"))
```

### Save Trees

```python
# Write single tree
Phylo.write(tree, "output.nwk", "newick")

# Write multiple trees
Phylo.write(trees, "all_trees.nex", "nexus")
```

### Format Conversion

```python
Phylo.convert("input.nwk", "newick", "output.xml", "phyloxml")
```

## Tree Structure

### Tree Components

- **Clade**: Any node (internal or leaf)
- **Terminal**: Leaf node (taxon/species)
- **Non-terminal**: Internal node
- **Branch length**: Evolutionary distance

### Accessing Nodes

```python
# Root clade
root = tree.root

# All terminal nodes (leaves)
terminals = tree.get_terminals()
print(f"Taxa: {len(terminals)}")

# Internal nodes only
internals = tree.get_nonterminals()

# All clades
all_clades = list(tree.find_clades())
```

## Tree Traversal

### Iteration Modes

```python
# Pre-order (root first)
for clade in tree.find_clades(order="preorder"):
    print(clade.name)

# Post-order (leaves first)
for clade in tree.find_clades(order="postorder"):
    print(clade.name)

# Level-order (breadth-first)
for clade in tree.find_clades(order="level"):
    print(clade.name)
```

### Terminal Iteration

```python
for taxon in tree.get_terminals():
    print(f"{taxon.name}: branch length = {taxon.branch_length}")
```

### Search by Criteria

```python
# Find by name
target = tree.find_any(name="Homo_sapiens")

# Find by custom predicate
def long_branch(clade):
    return clade.branch_length and clade.branch_length > 0.1

long_branches = list(tree.find_clades(long_branch))
```

## Tree Statistics

```python
# Total branch length
total_length = tree.total_branch_length()

# Depth (distance from root to tips)
depths = tree.depths()
max_depth = max(depths.values())

# Number of taxa
n_taxa = tree.count_terminals()
```

## Distance Calculations

### Between Taxa

```python
dist = tree.distance("Species_A", "Species_B")
print(f"Distance: {dist:.4f}")
```

### Distance Matrix

```python
def build_distance_matrix(tree):
    """Generate pairwise distance matrix."""
    taxa = [t.name for t in tree.get_terminals()]
    n = len(taxa)
    matrix = [[0.0] * n for _ in range(n)]

    for i, t1 in enumerate(taxa):
        for j, t2 in enumerate(taxa):
            if i < j:
                d = tree.distance(t1, t2)
                matrix[i][j] = d
                matrix[j][i] = d

    return taxa, matrix
```

## Common Ancestors

```python
# Find MRCA of two taxa
clade1 = tree.find_any(name="Species_A")
clade2 = tree.find_any(name="Species_B")
mrca = tree.common_ancestor(clade1, clade2)

# MRCA of multiple taxa
targets = [tree.find_any(name=n) for n in ["Sp_A", "Sp_B", "Sp_C"]]
mrca = tree.common_ancestor(*targets)
```

## Tree Manipulation

### Pruning (Remove Taxa)

```python
# Create copy before modifying
tree_copy = tree.copy()

# Remove single taxon
tree_copy.prune("Outgroup_species")

# Keep only specific taxa
taxa_to_keep = {"Species_A", "Species_B", "Species_C"}
for terminal in tree_copy.get_terminals():
    if terminal.name not in taxa_to_keep:
        tree_copy.prune(terminal)
```

### Rooting

```python
# Midpoint rooting
tree.root_at_midpoint()

# Outgroup rooting
outgroup = tree.find_any(name="Outgroup_sp")
tree.root_with_outgroup(outgroup)
```

### Ladderizing

```python
# Sort branches by descending size (largest first)
tree.ladderize(reverse=True)

# Sort by ascending size
tree.ladderize()
```

### Collapse Short Branches

```python
def collapse_branches(tree, threshold=0.001):
    """Set branches shorter than threshold to zero."""
    for clade in tree.find_clades():
        if clade.branch_length and clade.branch_length < threshold:
            clade.branch_length = 0
```

## Tree Visualization

### ASCII Display

```python
Phylo.draw_ascii(tree)
```

### Matplotlib Rendering

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 8))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_title("Phylogenetic Tree")
fig.tight_layout()
fig.savefig("tree.png", dpi=200)
plt.show()
```

### Display Options

```python
# Show branch lengths
Phylo.draw(tree, branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length else "")

# Custom label function
Phylo.draw(tree, label_func=lambda c: c.name.replace("_", " ") if c.name else "")
```

## Building Trees

### From Distance Matrix

```python
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

# Create distance matrix
names = ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"]
dm = DistanceMatrix(
    names=names,
    matrix=[
        [],
        [0.1],
        [0.3, 0.25],
        [0.5, 0.45, 0.2]
    ]
)

constructor = DistanceTreeConstructor()

# UPGMA tree
upgma_tree = constructor.upgma(dm)

# Neighbor-joining tree
nj_tree = constructor.nj(dm)
```

### From Sequence Alignment

```python
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Load alignment
alignment = AlignIO.read("aligned_seqs.fasta", "fasta")

# Calculate distances
calc = DistanceCalculator("identity")  # or "blosum62" for proteins
distance_matrix = calc.get_distance(alignment)

# Build tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

# Optional: root at midpoint
tree.root_at_midpoint()

Phylo.write(tree, "output.nwk", "newick")
```

### Distance Models

| Model | Use |
|-------|-----|
| `identity` | Simple identity distance |
| `blastn` | BLASTN-style scoring |
| `trans` | Transition/transversion |
| `blosum62` | BLOSUM62 for proteins |
| `pam250` | PAM250 for proteins |

## Consensus Trees

```python
from Bio.Phylo.Consensus import majority_consensus, strict_consensus

# Load bootstrap replicates
trees = list(Phylo.parse("bootstrap.nwk", "newick"))

# Majority-rule consensus (50% cutoff)
consensus = majority_consensus(trees, cutoff=0.5)

# Strict consensus (100% agreement)
strict = strict_consensus(trees)

Phylo.write(consensus, "consensus.nwk", "newick")
```

## Bootstrap Support

### Add Support Values

```python
def add_bootstrap_support(tree, support_dict):
    """Add bootstrap values to internal nodes."""
    for clade in tree.get_nonterminals():
        if clade.name in support_dict:
            clade.confidence = support_dict[clade.name]
```

### Display Support

```python
# Show confidence values
Phylo.draw(tree, label_func=lambda c: f"{c.confidence:.0f}" if c.confidence else "")
```

## Analysis Functions

### Extract Subtree

```python
def extract_clade(tree, taxa_names):
    """Extract subtree containing only specified taxa."""
    subtree = tree.copy()
    keep = set(taxa_names)

    for terminal in subtree.get_terminals():
        if terminal.name not in keep:
            subtree.prune(terminal)

    return subtree
```

### Phylogenetic Diversity

```python
def phylogenetic_diversity(tree, taxa=None):
    """Sum of branch lengths connecting specified taxa."""
    if taxa:
        tree = extract_clade(tree.copy(), taxa)

    total = sum(
        c.branch_length for c in tree.find_clades()
        if c.branch_length
    )
    return total
```

### Robinson-Foulds Distance

```python
def rf_distance(tree1, tree2):
    """Calculate Robinson-Foulds topological distance."""
    def get_splits(tree):
        splits = set()
        for clade in tree.get_nonterminals():
            leaves = frozenset(t.name for t in clade.get_terminals())
            splits.add(leaves)
        return splits

    splits1 = get_splits(tree1)
    splits2 = get_splits(tree2)

    return len(splits1.symmetric_difference(splits2))
```

### Compare Topologies

```python
def same_topology(tree1, tree2):
    """Check if two trees have identical topology."""
    # Same taxa
    taxa1 = set(t.name for t in tree1.get_terminals())
    taxa2 = set(t.name for t in tree2.get_terminals())

    if taxa1 != taxa2:
        return False

    return rf_distance(tree1, tree2) == 0
```

## PhyloXML Features

```python
from Bio.Phylo.PhyloXML import Phylogeny, Clade, Taxonomy

# Create annotated tree
tree = Phylogeny(rooted=True, name="Species Tree")

# Add taxonomy
clade = Clade(name="Homo_sapiens", branch_length=0.05)
taxonomy = Taxonomy(
    scientific_name="Homo sapiens",
    common_name="Human",
    rank="species"
)
clade.taxonomies.append(taxonomy)
```

## Recommendations

1. **Use Newick** for simple trees; phyloXML for annotations
2. **Validate trees** - check for polytomies, negative lengths
3. **Root appropriately** - midpoint or known outgroup
4. **Use tree.copy()** before modifications
5. **Document tree construction** - alignment, model, parameters
6. **Handle multiple trees** (bootstraps) with consensus methods
7. **Export publication-quality** figures with matplotlib
8. **Ensure consistent taxon naming** across files
9. **Consider tree size** - large trees may need special handling
10. **Validate bootstrap support** thresholds for confidence
