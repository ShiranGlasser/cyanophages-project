# from record.constants import CODON_TABLE, AMINO_TABLE

def calculate(dna_sequence, codon_table, amino_table):
        #assumes that dna_sequence is valid and a multiple of 3 

    codon_count = {}
    for k in range(0, len(dna_sequence), 3):
        codon = dna_sequence[k:k+3]
        amino = codon_table[codon]
        if amino in codon_count:
            if codon in codon_count[amino]:
                codon_count[amino][codon] += 1
            else:
                codon_count[amino][codon] = 1
        else:
            codon_count[amino] = {codon: 1}
    print(codon_count)
    
    Nk = 0
    
    for k in range(1, 7):
        if k==5:
            continue

        aminos_i = [f for f in amino_table.keys() if amino_table[f] == k and f != 'END']

        for amino_i in aminos_i:
            sigma_pi = 0
            n=0
                 
            if amino_i not in codon_count.keys():
                Fi = 1/k
                Nk += 1/Fi
                continue
          

            amino_i_dic = codon_count[amino_i]
            n  +=sum(amino_i_dic.values())
            if n==1:
                Fi = 1
                Nk += Fi
                continue
          
            codons_num = len(amino_i_dic)
            
            amino_i_codon_count = [val for val in amino_i_dic.values()]
            for j in range(0, codons_num):
                # import pdb;pdb.set_trace()
                pi =  amino_i_codon_count[j]/n
                sigma_pi += pow(pi,2)
            Fi = (n*sigma_pi-1)/(n-1)
            if Fi==0:
                Nk=k
            Nk+=1/Fi 

     
    return(Nk)
        

def stats_amino_acids(codon_table):
    dic = {}
    for elemnt in codon_table:
        if codon_table[elemnt] in dic.keys():
            dic[codon_table[elemnt]] = dic[codon_table[elemnt]] + 1
        elif codon_table[elemnt] == '_':
            if 'END' in dic.keys():
                dic['END'] =  dic['END'] + 1
            else:
                dic['END'] = 1
        else:
            dic[codon_table[elemnt]] = 1
    return dic

# dna_sequence = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGT"
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
AMINO_TABLE = {'I': 3, 'M': 1, 'T': 4, 'N': 2, 'K': 2, 'S': 6, 'R': 6, 'L': 6, 'P': 4, 'H': 2,
               'Q': 2, 'V': 4, 'A': 4, 'D': 2, 'E': 2, 'G': 4, 'F': 2, 'Y': 2, 'END': 3, 'C': 2, 'W': 1}
s='ATAATGACAAACAAAAGTAGACTCCCGCACCAACGGGTAGCAGATGAAGGGTTTTACTGCTGG' #20
s2="AAGAAA"
# print(calculate(s+s+s+s+s+s,codon_table,AMINO_TABLE))
sequence2 ="ATAATGACAAACAAAAGTAGACTGCCGCATCAAGTGGCTGATGAAGGGTTTTATTGTTGGATAATGACAAACAAAAGTAGACTGCCGCATCAAGTGGCTGATGAAGGGTTTTATTGTTGG"
sequence = "ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"
S="ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"
print(len(sequence)/3)
print((calculate(S,codon_table,AMINO_TABLE)))
# ATGGCCATGGCCACTGAGATCAATAGTACTGAGATCAATAGTGCT
# ATG TGC TCC CAG CCA TCA AGA AGT GTG AGT 
# M   C    S   Q   P   S  R   S   V    S   

# ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAG
#     '':'P', '':'P', '':'P', '':'P',
#     '':'H', '':'H', '':'Q', '':'Q',
#     'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
#     'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
#     'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
#     'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
#     'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
#     'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
#     'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
#     'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
#     'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',