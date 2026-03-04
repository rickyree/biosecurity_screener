# UniProt API Reference

Comprehensive technical reference for programmatic access to the UniProt protein database.

## Available APIs

| API | Purpose | Base URL |
|-----|---------|----------|
| Search API | Find entries matching criteria | `https://rest.uniprot.org/uniprotkb/search` |
| Data API | Retrieve entry by accession | `https://rest.uniprot.org/uniprotkb/{accession}` |
| ID Mapping | Translate identifiers | `https://rest.uniprot.org/idmapping/` |
| Stream API | Download large datasets | `https://rest.uniprot.org/uniprotkb/stream` |

## Query Syntax

### Boolean Operators

| Operator | Example | Description |
|----------|---------|-------------|
| `AND` | `kinase AND human` | Both conditions required |
| `OR` | `insulin OR glucagon` | Either condition |
| `NOT` | `kinase NOT mouse` | Exclude matches |
| `()` | `(A OR B) AND C` | Control precedence |

### Field Prefixes

#### Identity Fields
| Field | Example |
|-------|---------|
| `accession` | `accession:P00533` |
| `id` | `id:EGFR_HUMAN` |
| `gene` | `gene:EGFR` |
| `gene_exact` | `gene_exact:EGFR` |

#### Taxonomy Fields
| Field | Example |
|-------|---------|
| `organism_name` | `organism_name:"Homo sapiens"` |
| `organism_id` | `organism_id:9606` |
| `taxonomy_id` | `taxonomy_id:10090` |
| `taxonomy_name` | `taxonomy_name:"Mammalia"` |
| `lineage` | `lineage:vertebrata` |

#### Protein Information
| Field | Example |
|-------|---------|
| `protein_name` | `protein_name:kinase` |
| `recommended_name` | `recommended_name:insulin` |
| `reviewed` | `reviewed:true` |
| `fragment` | `fragment:false` |

#### Sequence Properties
| Field | Example |
|-------|---------|
| `length` | `length:[100 TO 500]` |
| `mass` | `mass:[20000 TO 50000]` |
| `sequence` | `sequence:MKTAYIAKQRQIS` |

#### Annotations
| Field | Example |
|-------|---------|
| `cc_function` | `cc_function:*` |
| `cc_disease` | `cc_disease:cancer` |
| `cc_subcellular_location` | `cc_subcellular_location:membrane` |
| `ft_signal` | `ft_signal:*` |
| `ft_transmem` | `ft_transmem:*` |
| `ft_domain` | `ft_domain:"SH2"` |

#### Cross-References
| Field | Example |
|-------|---------|
| `xref:pdb` | Has PDB structure |
| `xref:ensembl` | Has Ensembl link |
| `database:drugbank` | Has DrugBank entry |

#### Gene Ontology
| Field | Example |
|-------|---------|
| `go` | `go:0016301` (kinase activity) |
| `go_f` | `go_f:*` (molecular function) |
| `go_p` | `go_p:*` (biological process) |
| `go_c` | `go_c:*` (cellular component) |

### Range Expressions

```
length:[100 TO 500]       # 100-500 residues (inclusive)
mass:[* TO 30000]         # up to 30 kDa
mass:[50000 TO *]         # over 50 kDa
created:[2024-01-01 TO *] # created after Jan 2024
```

### Wildcards

```
gene:IL*              # IL1, IL2, IL6, IL10, etc.
protein_name:transport*
gene:BRCA?            # BRCA1, BRCA2
```

### Complex Query Examples

```
# Well-characterized human kinases with structure
(protein_name:kinase OR family:kinase) AND organism_id:9606 AND reviewed:true AND xref:pdb

# Secreted proteins with signal peptides
ft_signal:* AND cc_subcellular_location:secreted AND reviewed:true

# Disease proteins excluding mice
cc_disease:* AND organism_id:9606 NOT organism_name:mouse

# Recently updated entries
organism_id:9606 AND modified:[2024-01-01 TO *] AND reviewed:true

# Enzymes with known structure
cc_catalytic_activity:* AND xref:pdb AND reviewed:true
```

---

## API Fields Reference

Specify which data columns to retrieve using the `fields` parameter. Multiple fields are comma-separated.

```
https://rest.uniprot.org/uniprotkb/search?query=hemoglobin&fields=accession,gene_names,length
```

### Identification Fields

| Field | Description |
|-------|-------------|
| `accession` | Primary accession (e.g., P04637) |
| `id` | Entry name (e.g., P53_HUMAN) |
| `uniprotkb_id` | Same as id |
| `entryType` | REVIEWED (Swiss-Prot) or UNREVIEWED (TrEMBL) |

### Protein and Gene Names

| Field | Description |
|-------|-------------|
| `protein_name` | Recommended and alternative names |
| `gene_names` | All gene names |
| `gene_primary` | Primary gene symbol |
| `gene_synonym` | Gene name synonyms |
| `gene_oln` | Ordered locus names |
| `gene_orf` | ORF names |

### Organism Data

| Field | Description |
|-------|-------------|
| `organism_name` | Scientific name |
| `organism_id` | NCBI taxonomy ID |
| `lineage` | Full taxonomic lineage |
| `virus_hosts` | Host organisms (viral proteins) |

### Sequence Properties

| Field | Description |
|-------|-------------|
| `sequence` | Amino acid sequence |
| `length` | Residue count |
| `mass` | Molecular weight (Da) |
| `fragment` | Fragment indicator |
| `checksum` | CRC64 checksum |

### Functional Annotations (cc_)

#### Biochemistry
| Field | Description |
|-------|-------------|
| `cc_function` | Functional description |
| `cc_catalytic_activity` | Enzyme activity |
| `cc_activity_regulation` | Regulatory mechanisms |
| `cc_cofactor` | Required cofactors |
| `cc_pathway` | Metabolic pathway involvement |

#### Localization and Interactions
| Field | Description |
|-------|-------------|
| `cc_subcellular_location` | Cellular compartment |
| `cc_interaction` | Protein binding partners |
| `cc_subunit` | Quaternary structure |
| `cc_tissue_specificity` | Expression patterns |
| `cc_developmental_stage` | Developmental expression |

#### Medical Relevance
| Field | Description |
|-------|-------------|
| `cc_disease` | Associated diseases |
| `cc_disruption_phenotype` | Knockout phenotypes |
| `cc_allergen` | Allergenicity data |
| `cc_toxic_dose` | Toxicity information |

#### Modifications
| Field | Description |
|-------|-------------|
| `cc_ptm` | Post-translational modifications |
| `cc_mass_spectrometry` | MS evidence |
| `cc_alternative_products` | Isoforms/splice variants |
| `cc_polymorphism` | Polymorphism data |

### Sequence Features (ft_)

#### Processing Signals
| Field | Description |
|-------|-------------|
| `ft_signal` | Signal peptide |
| `ft_transit` | Transit peptide |
| `ft_propep` | Propeptide |
| `ft_chain` | Mature polypeptide |
| `ft_peptide` | Released peptide |

#### Structural Elements
| Field | Description |
|-------|-------------|
| `ft_domain` | Protein domain |
| `ft_repeat` | Repeated sequence |
| `ft_coiled` | Coiled-coil region |
| `ft_motif` | Short motif |
| `ft_region` | Region of interest |

#### Binding and Catalysis
| Field | Description |
|-------|-------------|
| `ft_act_site` | Catalytic residue |
| `ft_binding` | Ligand binding site |
| `ft_metal` | Metal coordination |
| `ft_dna_bind` | DNA interaction |
| `ft_np_bind` | Nucleotide binding |

#### Modifications
| Field | Description |
|-------|-------------|
| `ft_mod_res` | Modified residue |
| `ft_lipid` | Lipid attachment |
| `ft_carbohyd` | Glycosylation site |
| `ft_disulfid` | Disulfide bond |

#### Membrane Topology
| Field | Description |
|-------|-------------|
| `ft_transmem` | Transmembrane segment |
| `ft_intramem` | Intramembrane region |
| `ft_topo_dom` | Topological domain |

#### Variation
| Field | Description |
|-------|-------------|
| `ft_variant` | Natural variant |
| `ft_var_seq` | Alternative sequence |
| `ft_mutagen` | Mutagenesis site |
| `ft_conflict` | Sequence conflict |

### Gene Ontology

| Field | Description |
|-------|-------------|
| `go` | All GO annotations |
| `go_p` | Biological process |
| `go_c` | Cellular component |
| `go_f` | Molecular function |
| `go_id` | GO identifiers |

### Database Cross-References (xref_)

#### Sequence Resources
| Field | Database |
|-------|----------|
| `xref_embl` | EMBL/GenBank/DDBJ |
| `xref_refseq` | RefSeq |
| `xref_ccds` | CCDS |

#### Structure Resources
| Field | Database |
|-------|----------|
| `xref_pdb` | Protein Data Bank |
| `xref_alphafolddb` | AlphaFold DB |
| `xref_smr` | SWISS-MODEL |

#### Protein Families
| Field | Database |
|-------|----------|
| `xref_interpro` | InterPro |
| `xref_pfam` | Pfam |
| `xref_prosite` | PROSITE |
| `xref_smart` | SMART |

#### Genomic Resources
| Field | Database |
|-------|----------|
| `xref_ensembl` | Ensembl |
| `xref_geneid` | NCBI Gene |
| `xref_kegg` | KEGG |

#### Disease Resources
| Field | Database |
|-------|----------|
| `xref_omim` | OMIM |
| `xref_orphanet` | Orphanet |
| `xref_disgenet` | DisGeNET |

#### Pharmacology
| Field | Database |
|-------|----------|
| `xref_chembl` | ChEMBL |
| `xref_drugbank` | DrugBank |
| `xref_guidetopharmacology` | GtoPdb |

### Entry Metadata

| Field | Description |
|-------|-------------|
| `date_created` | Creation date |
| `date_modified` | Last update |
| `date_sequence_modified` | Sequence update |
| `annotation_score` | Quality score (1-5) |
| `protein_existence` | Evidence level |
| `reviewed` | Swiss-Prot status |
| `lit_pubmed_id` | PubMed references |

### Common Field Sets

**Minimal identification:**
```
accession,id,protein_name,gene_names,organism_name
```

**Sequence analysis:**
```
accession,sequence,length,mass,xref_pdb,xref_alphafolddb
```

**Functional profiling:**
```
accession,protein_name,cc_function,cc_catalytic_activity,go,cc_pathway
```

**Clinical applications:**
```
accession,gene_names,cc_disease,xref_omim,xref_malacards,ft_variant
```

**Structural biology:**
```
accession,sequence,xref_pdb,xref_alphafolddb,ft_domain,ft_transmem
```

### Field Discovery

```python
import requests
resp = requests.get("https://rest.uniprot.org/configure/uniprotkb/result-fields")
all_fields = resp.json()
```

---

## ID Mapping Databases

Supported database names for the UniProt ID Mapping service. Use exact names in API calls.

### UniProt Resources

| Database | Name |
|----------|------|
| UniProt accession/ID | `UniProtKB_AC-ID` |
| UniProtKB (all) | `UniProtKB` |
| Swiss-Prot (reviewed) | `UniProtKB-Swiss-Prot` |
| TrEMBL (unreviewed) | `UniProtKB-TrEMBL` |
| UniProt Archive | `UniParc` |
| UniRef 50% clusters | `UniRef50` |
| UniRef 90% clusters | `UniRef90` |
| UniRef 100% clusters | `UniRef100` |

### Nucleotide Sequences

| Database | Name |
|----------|------|
| EMBL/GenBank/DDBJ | `EMBL` |
| EMBL CDS | `EMBL-CDS` |
| RefSeq nucleotide | `RefSeq_Nucleotide` |
| CCDS | `CCDS` |

### Protein Sequences

| Database | Name |
|----------|------|
| RefSeq protein | `RefSeq_Protein` |
| PIR | `PIR` |

### Gene Identifiers

| Database | Name |
|----------|------|
| NCBI Gene | `GeneID` |
| Gene name | `Gene_Name` |
| Gene synonym | `Gene_Synonym` |
| Ordered locus name | `Gene_OrderedLocusName` |
| ORF name | `Gene_ORFName` |

### Genome Browsers

| Database | Name |
|----------|------|
| Ensembl gene | `Ensembl` |
| Ensembl protein | `Ensembl_PRO` |
| Ensembl transcript | `Ensembl_TRS` |
| KEGG Genes | `KEGG` |
| UCSC | `UCSC` |

### 3D Structure

| Database | Name |
|----------|------|
| PDB | `PDB` |
| AlphaFold DB | `AlphaFoldDB` |
| SWISS-MODEL | `SMR` |

### Protein Families and Domains

| Database | Name |
|----------|------|
| InterPro | `InterPro` |
| Pfam | `Pfam` |
| PROSITE | `PROSITE` |
| SMART | `SMART` |
| PANTHER | `PANTHER` |

### Model Organisms

| Database | Organism | Name |
|----------|----------|------|
| MGI | Mouse | `MGI` |
| RGD | Rat | `RGD` |
| FlyBase | Drosophila | `FlyBase` |
| WormBase | C. elegans | `WormBase` |
| SGD | Yeast (S. cerevisiae) | `SGD` |
| TAIR | Arabidopsis | `TAIR` |
| HGNC | Human | `HGNC` |

### Pathways and Metabolism

| Database | Name |
|----------|------|
| Reactome | `Reactome` |
| BioCyc | `BioCyc` |
| SIGNOR | `SIGNOR` |

### Disease and Phenotype

| Database | Name |
|----------|------|
| OMIM | `OMIM` |
| Orphanet | `OrphaNet` |
| DisGeNET | `DisGeNET` |
| Open Targets | `OpenTargets` |

### Pharmacology

| Database | Name |
|----------|------|
| ChEMBL | `ChEMBL` |
| DrugBank | `DrugBank` |
| DrugCentral | `DrugCentral` |

### Protein Interactions

| Database | Name |
|----------|------|
| STRING | `STRING` |
| BioGRID | `BioGRID` |
| IntAct | `IntAct` |

### Orthology and Evolution

| Database | Name |
|----------|------|
| Gene Ontology | `GO` |
| OMA | `OMA` |
| OrthoDB | `OrthoDB` |
| eggNOG | `eggNOG` |

### Discover Available Databases

```
GET https://rest.uniprot.org/configure/idmapping/fields
```

---

## Code Examples

### curl

**Basic search:**
```bash
curl "https://rest.uniprot.org/uniprotkb/search?query=hemoglobin&format=json&size=3"
```

**Fetch single entry:**
```bash
curl "https://rest.uniprot.org/uniprotkb/P04637.fasta"
```

**Custom fields (TSV):**
```bash
curl "https://rest.uniprot.org/uniprotkb/search?query=gene:TP53&format=tsv&fields=accession,gene_names,length"
```

**Submit ID mapping:**
```bash
curl -X POST "https://rest.uniprot.org/idmapping/run" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "from=UniProtKB_AC-ID&to=PDB&ids=P04637,P00533"
```

**Download complete dataset:**
```bash
curl "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta" \
  -o human_reviewed.fasta
```

### R

**Search and parse JSON:**
```r
library(httr)
library(jsonlite)

endpoint <- "https://rest.uniprot.org/uniprotkb/search"
params <- list(
  query = "hemoglobin AND organism_id:9606",
  format = "json",
  size = 5
)

resp <- GET(endpoint, query = params)
data <- fromJSON(content(resp, "text"))

proteins <- data.frame(
  accession = data$results$primaryAccession,
  stringsAsFactors = FALSE
)
print(proteins)
```

**TSV to data frame:**
```r
library(httr)
library(readr)

params <- list(
  query = "gene:EGFR AND reviewed:true",
  format = "tsv",
  fields = "accession,gene_names,organism_name,length"
)

resp <- GET("https://rest.uniprot.org/uniprotkb/search", query = params)
df <- read_tsv(content(resp, "text"))
print(df)
```

### JavaScript

**Fetch API search:**
```javascript
async function searchUniProt(query) {
  const params = new URLSearchParams({
    query: query,
    format: "json",
    size: "10"
  });

  const resp = await fetch(
    `https://rest.uniprot.org/uniprotkb/search?${params}`
  );
  return resp.json();
}

const data = await searchUniProt("hemoglobin AND organism_id:9606");
console.log(data.results);
```

**ID mapping workflow:**
```javascript
async function mapIdentifiers(ids, fromDb, toDb) {
  // Submit
  const submitResp = await fetch("https://rest.uniprot.org/idmapping/run", {
    method: "POST",
    body: new URLSearchParams({
      from: fromDb,
      to: toDb,
      ids: ids.join(",")
    })
  });
  const { jobId } = await submitResp.json();

  // Poll
  while (true) {
    const statusResp = await fetch(
      `https://rest.uniprot.org/idmapping/status/${jobId}`
    );
    const status = await statusResp.json();
    if ("results" in status || "failedIds" in status) break;
    await new Promise(r => setTimeout(r, 2000));
  }

  // Fetch results
  const resultsResp = await fetch(
    `https://rest.uniprot.org/idmapping/results/${jobId}`
  );
  return resultsResp.json();
}
```

---

## Rate Limiting

UniProt enforces rate limits to ensure fair access.

### Guidelines

- Start with 2-3 requests per second
- Implement exponential backoff on 429 responses
- Cache results to avoid redundant requests
- Use streaming for large result sets

### Backoff Implementation

```python
import time
import requests

def fetch_with_backoff(url, max_retries=5):
    """Fetch URL with automatic retry on rate limits."""
    delay = 1.0

    for attempt in range(max_retries):
        resp = requests.get(url)

        if resp.status_code == 200:
            return resp
        elif resp.status_code == 429:
            print(f"Rate limited, waiting {delay}s...")
            time.sleep(delay)
            delay *= 2
        else:
            resp.raise_for_status()

    raise RuntimeError(f"Failed after {max_retries} attempts")
```

---

## Error Handling

| HTTP Code | Meaning | Action |
|-----------|---------|--------|
| 200 | Success | Process response |
| 400 | Bad request | Check query syntax |
| 404 | Entry not found | Verify accession, may be obsolete |
| 429 | Rate limited | Wait and retry with backoff |
| 500 | Server error | Retry after brief delay |

### Robust Request Function

```python
import requests
import time

def safe_request(url, params=None, max_retries=3):
    """Make request with error handling and retries."""
    delay = 1.0

    for attempt in range(max_retries):
        try:
            resp = requests.get(url, params=params, timeout=30)

            if resp.status_code == 200:
                return resp
            elif resp.status_code == 429:
                time.sleep(delay)
                delay *= 2
            elif resp.status_code == 404:
                return None  # Entry doesn't exist
            else:
                resp.raise_for_status()

        except requests.RequestException as e:
            if attempt == max_retries - 1:
                raise
            time.sleep(delay)
            delay *= 2

    return None
```

---

## Notes

- Database names in ID mapping are case-sensitive
- One-to-many mappings are common; check `failedIds` for unmapped identifiers
- Wildcards work for field groups: `cc_*` retrieves all comment fields
- Fewer fields = faster response times
- JSON returns structured objects; TSV returns flattened strings
- ID mapping results expire after 7 days
- Maximum 100,000 identifiers per mapping request

## Additional Resources

- REST API: https://www.uniprot.org/help/api
- Query Fields: https://www.uniprot.org/help/query-fields
- ID Mapping: https://www.uniprot.org/help/id_mapping
- SPARQL Endpoint: https://sparql.uniprot.org/
- API Interactive Docs: https://www.uniprot.org/api-documentation
