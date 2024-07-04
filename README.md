# Genetic Code Translator

This Python code translates a DNA sequence into its corresponding protein sequence.

### Functionality

* **DNA Codon Table:** Utilizes a dictionary (`codon_table`) to map DNA codons (triplets of nucleotides) to their corresponding amino acids.
* **Iterative Translation:** Iterates through the DNA sequence in steps of 3 (codon length) to extract codons.
* **Amino Acid Lookup:** Looks up each codon in the `codon_table` to retrieve the associated amino acid.
* **Stop Codon Handling:** Terminates translation upon encountering a stop codon ("*").
* **Unknown Codon Handling:** Represents unknown codons (not found in the table) with "X" in the protein sequence.
* **Protein Sequence Construction:** Builds the protein sequence by concatenating the translated amino acids.

### Example Usage

```python
dna_sequence = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGATGCGGCGCAGGAAGGGGTCGGAGTGA"
protein_sequence = translate_dna_to_protein(dna_sequence)
print("DNA Sequence:", dna_sequence)
print("Protein Sequence:", protein_sequence)
