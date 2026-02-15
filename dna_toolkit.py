# dna_toolkit.py
import collections

class DNA_Toolkit:
    def __init__(self, sequence):
        # This cleans the sequence: makes it uppercase and removes spaces
        self.sequence = sequence.upper().replace(" ", "")
        
        # Check if the sequence is valid immediately
        if not self._validate():
            raise ValueError("Invalid DNA Sequence: Contains characters other than A, T, G, C")

    def _validate(self):
        """Check if the sequence only contains A, T, G, C"""
        valid_bases = set("ATGC")
        # returns True if all characters in sequence are in valid_bases
        return set(self.sequence).issubset(valid_bases)

    def count_nucleotides(self):
        """Count the number of A, T, G, and C"""
        # Counter creates a dictionary like {'A': 10, 'T': 12...}
        return dict(collections.Counter(self.sequence))

    def gc_content(self):
        """Calculate GC Content percentage"""
        c_count = self.sequence.count('C')
        g_count = self.sequence.count('G')
        return (c_count + g_count) / len(self.sequence) * 100

    def transcribe(self):
        """Transcribes DNA to RNA (Replaces T with U)"""
        return self.sequence.replace('T', 'U')
    
    def reverse_complement(self):
        """Generates the reverse complement of the DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        # Create the reverse complement by mapping each base and reversing the string
        return ''.join(complement[base] for base in reversed(self.sequence))
    
    def is_palindrome(self):
        """Checks if the DNA sequence is a palindrome (same as its reverse complement)"""
        return self.sequence == self.reverse_complement()
    
    def find_enzyme_sites(self, enzyme_seq="GAATTC"):
        """Finds index positions of a restriction enzyme (Default: EcoRI)"""
        sites = []
        index = self.sequence.find(enzyme_seq)
        
        while index != -1:
            sites.append(index)
            # Search again starting from the next character
            index = self.sequence.find(enzyme_seq, index + 1)
            
        return sites if sites else "No sites found"
    
    def protein_weight(self):
        """Calculates total molecular weight of the protein in Daltons (Da)"""
        protein_seq = self.translate()
        
        # Weights of Amino Acids (average isotopic masses)
        aa_weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
            'Q': 146.15, 'E': 147.13, 'G': 75.07,  'H': 155.16, 'I': 131.18,
            'L': 131.18, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
            '_': 0.0, 'X': 0.0 # Stop codon & unknown
        }
        
        total_weight = sum(aa_weights.get(aa, 0) for aa in protein_seq)
        return round(total_weight, 2)

    def translate(self):
        """Translates RNA to Protein using a Codon Table"""
        # 1. Get RNA sequence first
        rna = self.transcribe()
        
        # 2. Define the Codon Table (The Genetic Code)
        codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }

        protein = ""
        # 3. Read sequence in chunks of 3 (Codons)
        for i in range(0, len(rna) - 2, 3):
            codon = self.sequence[i:i+3] # Use DNA codon table logic directly
            amino_acid = codon_table.get(codon, 'X') # 'X' if unknown
            
            if amino_acid == '_': # Stop codon
                break
            protein += amino_acid
            
        return protein

# --- MAIN EXECUTION BLOCK ---
if __name__ == "__main__":
    # A longer sequence that includes an EcoRI site (GAATTC)
    test_dna = "agccctccag gacaggctgc atcagaagag gccatcaagc agatcactgt ccttctgcca  tggccctgtg gatgcgcctc ctgcccctgc tggcgctgct ggccctctgg ggacctgacc cagccgcagc ctttgtgaac caacacctgt gcggctcaca cctggtggaa gctctctacctagtgtgcgg ggaacgaggc ttcttctaca cacccaagac ccgccgggag gcagaggacctgcaggtggg gcaggtggag ctgggcgggg gccctggtgc aggcagcctg cagcccttggccctggaggg gtccctgcag aagcgtggca ttgtggaaca atgctgtacc agcatctgctccctctacca gctggagaac tactgcaact agacgcagcc cgcaggcagc cccacacccgcgcctcctg caccgagaga gatggaataa agcccttgaa ccagcTAG" 
    
    print(f"--- DNA Toolkit Analysis ---")
    print(f"Input Sequence: {test_dna}\n")

    tool = DNA_Toolkit(test_dna)

    print(f"1. Nucleotide Count:  {tool.count_nucleotides()}")
    print(f"2. GC Content:        {tool.gc_content():.2f}%")
    print(f"3. Reverse Comp:      {tool.reverse_complement()}")
    print(f"8. reverse_complement {tool.reverse_complement()}")
    print(f"4. Is Palindrome?     {'Yes' if tool.is_palindrome() else 'No'}")
    print(f"5. Protein Sequence:  {tool.translate()}")
    print(f"6. Protein Weight:    {tool.protein_weight()} Da")
    print(f"7. EcoRI Cut Sites:   Found at index {tool.find_enzyme_sites('GAATTC')}")