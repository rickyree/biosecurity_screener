import random

codon_table = {
    'A': ['GCT','GCC','GCA','GCG'],
    'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
    'N': ['AAT','AAC'],
    'D': ['GAT','GAC'],
    'C': ['TGT','TGC'],
    'Q': ['CAA','CAG'],
    'E': ['GAA','GAG'],
    'G': ['GGT','GGC','GGA','GGG'],
    'H': ['CAT','CAC'],
    'I': ['ATT','ATC','ATA'],
    'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
    'K': ['AAA','AAG'],
    'M': ['ATG'],
    'F': ['TTT','TTC'],
    'P': ['CCT','CCC','CCA','CCG'],
    'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
    'T': ['ACT','ACC','ACA','ACG'],
    'W': ['TGG'],
    'Y': ['TAT','TAC'],
    'V': ['GTT','GTC','GTA','GTG']
}

input_tsv = "candidate_sequences.tsv"
output_fasta = "cholix_t1.0.fasta"

with open(input_tsv) as f:
    lines = f.readlines()

# rows 49–73 (inclusive) when header exists
selected = lines[302:361]

with open(output_fasta, "w") as out:
    for line in selected:
        cols = line.strip().split("\t")

        name = cols[0]
        protein_seq = cols[1]

        valid_aas = list(codon_table.keys())

        dna_seq = "".join(
            random.choice(codon_table[random.choice(valid_aas)] if aa == "X" else codon_table[aa])
            for aa in protein_seq
        )

        out.write(f">{name}\n")

        for i in range(0, len(dna_seq), 60):
            out.write(dna_seq[i:i+60] + "\n")