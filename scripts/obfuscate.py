import random

def process_dna(sequence):
    seq_len = len(sequence)
    min_len = 20
    
    # --- (b) Split into fragments ---
    # Determine how many fragments are possible (max 5)
    max_possible_frags = seq_len // min_len
    num_frags = min(5, max_possible_frags)
    
    if num_frags < 1:
        return "Sequence too short to process."

    # Generate random cut points
    # We need (num_frags - 1) unique points to create num_frags pieces
    cut_points = sorted(random.sample(range(min_len, seq_len - min_len + 1), num_frags - 1))
    
    # Ensure every fragment is at least 20 bases (Validation/Adjustment)
    # Adding 0 and the total length to the points to slice the string
    indices = [0] + cut_points + [seq_len]
    fragments = [sequence[indices[i]:indices[i+1]] for i in range(len(indices)-1)]

    processed_fragments = []
    
    # Helper for Reverse Complement
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    
    # --- (c) & (d) Assign sense/antisense and Scramble ---
    for frag in fragments:
        is_antisense = random.choice([True, False])
        
        if is_antisense:
            # Reverse complement the fragment
            rev_comp = "".join(complement.get(base, base) for base in reversed(frag))
            processed_fragments.append(rev_comp)
        else:
            processed_fragments.append(frag)
            
    # Scramble the order of the final fragments
    random.shuffle(processed_fragments)
    
    return "".join(processed_fragments)

# Example Usage:
target_prefix = "cholix"   # change this

sequences = {}

# file = "candidate_sequences_dna.fasta"
file = "cholix_t1.0.fasta"

with open(file) as f:
    header = None
    seq_lines = []

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            # save previous sequence
            if header and header.startswith(target_prefix):
                sequences[header] = "".join(seq_lines)

            header = line[1:]  # remove ">"
            seq_lines = []

        else:
            seq_lines.append(line)

    # save last sequence
    if header and header.startswith(target_prefix):
        sequences[header] = "".join(seq_lines)

import os

output_dir = "commec_seq/obfused/cholix_t1.0"

os.makedirs(output_dir, exist_ok=True)

for name, dna in sequences.items():
    result = process_dna(dna)

    # make filename safe
    safe_name = name.replace(" ", "_").replace("|", "_")

    filename = os.path.join(output_dir, f"{safe_name}.fasta")

    with open(filename, "w") as f:
        f.write(f">{safe_name}\n")
        f.write(result + "\n")