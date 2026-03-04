from Bio import Align

def perform_alignment(seq1, seq2, mode='global'):
    # 1. Initialize the aligner
    aligner = Align.PairwiseAligner()
    
    # 2. Set the alignment mode ('global' or 'local')
    aligner.mode = mode
    
    # 3. Configure scoring parameters (customize these as needed)
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    
    # 4. Perform the alignment
    alignments = aligner.align(seq1, seq2)
    
    # 5. Output the results
    print(f"--- {mode.capitalize()} Alignment ---")
    print(f"Number of optimal alignments: {len(alignments)}")
    print(f"Optimal Score: {alignments.score}\n")
    
    # Print the top alignment (or loop through 'alignments' to see all ties)
    if alignments:
        print("Top Alignment:")
        print(alignments[0])

# Example Sequences
sequence_A = "IGLSVINNGEYIVTTNATLYHGTKLGYTNSIQNGIQQPLNGTQGNYDTAHMGFYSTRPKFLAAGYAQDNINRFIQRRQGVAKVNFPGQTATRAHEIVNDEMRKIALLRTKDIPLMAIVGTSGFISVFGDRASSVVLTKSQAPGNSGMEHIMNWKKALYGMDKLELNQQKQEMPSGHAYYVAISRCCMCLNLNYEPIRKQAAQVMTNLMRYPAVSACTTENKQTVVWINNSKLRIAGLVEKILWHKSMINWKKMMSKEKNVSSQNQAELPSQCASLVNNAMSKNENKHMETLSTVGGWNTTFGYFNGAMYFNNEQILAQTRATCWYLVRQNMTLVGEKVNWASMEIIFTTSVIILDHINPERMNRVTQITEFQMVRIQADNFEIAMYSEVDRVILSRFRGSRIMRLLSSPETGNRPIRGILIPTLPELISVNNSQTQIAANQYATRFTCTHRTRMKIYCVPFNPVYLARPIHGGICIEVLRDSFEALNENTLASHKVVVLCGKLKWKNKWRFNYISITFVIRS"
sequence_B = "IGLSVINNGEYIVTTNATLYHGTKLGYTNSIQNGIQQPLNGTQGNYDTAHMGFYSTRPKFLAAGYAQDNINRFIQRRQGVAKVNFPGQTATRAHEIVNDEMRKIALLRTKDIPLMAIVGTSGFISVFGDRASSVVLTKSQAPGNSGMEHIMNWKKALYGMDKLELNQQKQEMPSGHAYYVAISRCCMXXXXXXXXXXXXXCLNLNYEPIRKQAAQVMTNLMRYPAVSACTTENKQTVVWINNSKLRIAGLVEKILWHKSMINWKKMMSKEKNVSSQNQAELPSQCASLVNNAMSKNENKHMETLSTVGGWNTTFGYFNGAMYFNNEQILAQTRATCWYLVRQNMTLVGEKVNWASMEIIFTTSVIILDHINPERMNRVTQITEFQMVRIQADNFEIAMYSEVDRVILSRFRGSRIMRLLSSPETGNRPIRGILIPTLPELISVNNSQTQIAANQYATRFTCTHRTRMKIYCVPFNPVYLARPIHGGICIEVLRDSFEALNENTLASHKVVVLCGKLKWKNKWRFNYISITFVIRS"

# Run Global Alignment (Needleman-Wunsch)
# Best for sequences of similar length that are expected to align end-to-end.
perform_alignment(sequence_A, sequence_B, mode='global')

print("\n" + "="*40 + "\n")

# Run Local Alignment (Smith-Waterman)
# Best for finding similar sub-regions between sequences of different lengths.
perform_alignment(sequence_A, sequence_B, mode='local') 
