# dna_analyzer.py

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence.upper()

def nucleotide_count(sequence):
    return {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }

def transcribe_dna_to_rna(sequence):
    return sequence.replace('T', 'U')

def gc_content(sequence):
    gc = sequence.count('G') + sequence.count('C')
    return (gc / len(sequence)) * 100

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([complement[base] for base in reversed(sequence)])

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

def translate_dna(sequence):
    protein = ''
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]
        protein += codon_table.get(codon, 'X')  # X = unknown codon
    return protein

if __name__ == "__main__":
    seq = read_fasta("sequences/example.fasta")
    print("Nucleotide Count:", nucleotide_count(seq))
    print("Transcribed RNA:", transcribe_dna_to_rna(seq))
    print("Translated Protein:", translate_dna(seq))
    print("GC Content: {:.2f}%".format(gc_content(seq)))
    print("Reverse Complement:", reverse_complement(seq))
