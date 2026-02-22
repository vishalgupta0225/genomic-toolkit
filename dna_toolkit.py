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
        return set(self.sequence).issubset(valid_bases)

    def count_nucleotides(self):
        """Count the number of A, T, G, and C"""
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
        return ''.join(complement[base] for base in reversed(self.sequence))
    
    def is_palindrome(self):
        """Checks if the DNA sequence is a palindrome"""
        return self.sequence == self.reverse_complement()
    
    def melting_temperature(self):
        """Calculates melting temperature (Tm)"""
        length = len(self.sequence)
        if length == 0: 
            return 0
        
        a_count = self.sequence.count('A')
        t_count = self.sequence.count('T')
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        
        if length < 14:
            tm = (2 * (a_count + t_count)) + (4 * (g_count + c_count))
        else:
            tm = 64.9 + 41 * (g_count + c_count - 16.4) / length
        return round(tm, 2)

    def find_all_restriction_sites(self):
        """Finds all common restriction sites in the sequence"""
        common_enzymes = {
            "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
            "NotI": "GCGGCCGC", "XhoI": "CTCGAG", "SmaI": "CCCGGG",
            "PstI": "CTGCAG", "TaqI": "TCGA"
        }
        results = {} 
        
        for enzyme_name, recognition_site in common_enzymes.items():
            positions = []
            for i in range(len(self.sequence) - len(recognition_site) + 1):
                if self.sequence[i:i+len(recognition_site)] == recognition_site:
                    positions.append(i + 1) # 1-based indexing for biology
            if positions:
                results[enzyme_name] = positions
                
        return results

    def protein_weight(self):
        """Calculates total molecular weight of the protein in Daltons (Da)"""
        protein_seq = self.translate()
        aa_weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
            'Q': 146.15, 'E': 147.13, 'G': 75.07,  'H': 155.16, 'I': 131.18,
            'L': 131.18, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
            '_': 0.0, 'X': 0.0 
        }
        total_weight = sum(aa_weights.get(aa, 0) for aa in protein_seq)
        return round(total_weight, 2)

    def translate(self):
        """Translates DNA to Protein using a Codon Table."""
        dna = self.sequence
        codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        protein = ""
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i+3]
            amino_acid = codon_table.get(codon, 'X')
            if amino_acid == '_':
                break
            protein += amino_acid
        return protein
    
    def evaluate_primer(self):
        """Evaluates a DNA sequence for use as a PCR primer."""
        length = len(self.sequence)
        gc_perc = self.gc_content()
        tm = self.melting_temperature()
        
        # Check if the very last letter is G or C (a GC Clamp)
        has_gc_clamp = False
        if length > 0:
            last_base = self.sequence[-1] # Gets the last character
            if last_base == 'G' or last_base == 'C':
                has_gc_clamp = True
                
        # Simple evaluation logic using standard biology rules
        results = {}
        
        # 1. Check Length (Ideal: 18 to 25 bp)
        if 18 <= length <= 25:
            results["Length"] = {"Value": f"{length} bp", "Status": "Good"}
        else:
            results["Length"] = {"Value": f"{length} bp", "Status": "Warning (Ideal is 18-25 bp)"}
            
        # 2. Check GC Content (Ideal: 40% to 60%)
        if 40 <= gc_perc <= 60:
            results["GC Content"] = {"Value": f"{gc_perc:.1f}%", "Status": "Good"}
        else:
            results["GC Content"] = {"Value": f"{gc_perc:.1f}%", "Status": "Warning (Ideal is 40-60%)"}
            
        # 3. Check Melting Temp (Ideal: 50°C to 65°C)
        if 50 <= tm <= 65:
            results["Melting Temp (Tm)"] = {"Value": f"{tm} °C", "Status": "Good"}
        else:
            results["Melting Temp (Tm)"] = {"Value": f"{tm} °C", "Status": "Warning (Ideal is 50-65°C)"}
            
        # 4. Check GC Clamp
        if has_gc_clamp:
            results["GC Clamp (Ends in G or C)"] = {"Value": "Yes", "Status": "Good"}
        else:
            results["GC Clamp (Ends in G or C)"] = {"Value": "No", "Status": "Warning (Should end in G or C)"}
            
        return results
