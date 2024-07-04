
def translate_dna_to_protein(dna_sequence):
    # Dictionary mapping DNA codons to amino acids
    codon_table = {
        "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",  # Isoleucine (I), Methionine (M)
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",  # Threonine (T)
        "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",  # Asparagine (N), Lysine (K)
        "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",  # Serine (S), Arginine (R)
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",  # Leucine (L)
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",  # Proline (P)
        "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",  # Histidine (H), Glutamine (Q)
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",  # Arginine (R)
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",  # Valine (V)
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",  # Alanine (A)
        "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",  # Aspartic Acid (D), Glutamic Acid (E)
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",  # Glycine (G)
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",  # Serine (S)
        "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",  # Phenylalanine (F), Leucine (L)
        "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",  # Tyrosine (Y), Stop (*)
        "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",   # Cysteine (C), Tryptophan (W), Stop (*)
    }

    protein_sequence = ""  # Initialize an empty string to store the protein sequence
    for i in range(0, len(dna_sequence), 3):  # Iterate over the DNA sequence with a step of 3 (codon length)
        codon = dna_sequence[i:i + 3]  # Extract a codon (sequence of 3 nucleotides)
        if codon in codon_table:  # Check if the codon exists in the codon table
            amino_acid = codon_table[codon]  # Get the corresponding amino acid for the codon
            if amino_acid == "*":  # Check if the amino acid is a stop codon
                break  # If it is a stop codon, terminate translation
            protein_sequence += amino_acid  # Add the amino acid to the protein sequence
        else:
            protein_sequence += "X"  # If the codon is not found in the table, represent it as 'X' (unknown codon)

    return protein_sequence  # Return the translated protein sequence

# Example usage:
dna_sequence = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGATGCGGCGCAGGAAGGGGTCGGAGTGA"
protein_sequence = translate_dna_to_protein(dna_sequence)
print("DNA Sequence:", dna_sequence)
print("Protein Sequence:", protein_sequence)




