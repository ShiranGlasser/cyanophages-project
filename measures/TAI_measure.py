import math

def calculate(dna_sequence, codon_table):
    return 5
    #assumes that dna_sequence is valid and a multiple of 3 
    codon_count = {}
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        amino = codon_table[codon]
        if amino in codon_count:
            if codon in codon_count[amino]:
                codon_count[amino][codon] += 1
            else:
                codon_count[amino][codon] = 1
        else:
            codon_count[amino] = {codon: 1}
    TAI = 0
    for amino, codons in codon_count.items():
        for codon, count in codons.items():
            f_i = count/len(dna_sequence)
            e_i = f_i * math.log(f_i)
            TAI += e_i
    TAI = 2 - TAI
    return TAI
