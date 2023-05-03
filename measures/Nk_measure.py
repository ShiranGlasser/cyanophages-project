from constants import AMINO_TABLE, CODON_TABLE

def calculate(dna_sequence):
    """
        This method calculates the Nk measure for the given dna_sequence
        The method assumes that dna_sequence is valid and a multiple of 3 
    """
    # possible future extensions written as TODO in the code

    # TODO 1- receive the pre calculated Faa for all the sequence's
    #           genome and use it when a certain amino in missing
    codon_count = {}
    for k in range(0, len(dna_sequence), 3):
        codon = dna_sequence[k:k+3]
        if len(codon)<3:
            break
        amino = CODON_TABLE[codon]
        if amino in codon_count:
            if codon in codon_count[amino]:
                codon_count[amino][codon] += 1
            else:
                codon_count[amino][codon] = 1
        else:
            codon_count[amino] = {codon: 1}
    # print(codon_count)

    Nk = 0

    for k in range(1, 7):
        if k == 5:
            continue

        aminos_i = [f for f in AMINO_TABLE.keys() if AMINO_TABLE[f]
                    == k and f != 'END']

        for amino_i in aminos_i:
            sigma_pi = 0
            n = 0

            if amino_i not in codon_count.keys():
                Fi = 1/k
                # TODO 2 - add option for Fi = the missing amino in the whole genome
                # TODO 3 - add option for Fi = the average of all the other aminos with the same type
                Nk += 1/Fi
                continue

            # containing for the current amino, its decoding codons and their appearances on the sequence
            amino_i_dic = codon_count[amino_i]
            # containing the number of deferent codons decoding the amino on the sequence
            codons_num = len(amino_i_dic)
            # the total number of the appearances of the current amino acid
            n += sum(amino_i_dic.values())
            synchronize="" 
            if n == 1 or codons_num == 1:
                Fi = 1
                Nk += Fi
                continue
            # the codons distribution for the amino is equal
            if len(set(amino_i_dic.values())) == 1:
                Fi = 1/k
                Nk += Fi
                continue

            amino_i_codon_count = [val for val in amino_i_dic.values()] # list of the numbers of performances for each codon 
            for j in range(0, codons_num):
                pi = amino_i_codon_count[j]/n
                sigma_pi += pow(pi, 2)
            Fi = (n*sigma_pi-1)/(n-1)
            if Fi == 0:
                print(
                    f"The case for this sequence:{dna_sequence}, is not yes supported")
                return -1
            Nk += 1/Fi

    return (Nk)



# dna_sequence = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGT"

s = 'ATAATGACAAACAAAAGTAGACTCCCGCACCAACGGGTAGCAGATGAAGGGTTTTACTGCTGG'  # 20
s2 = "AAGAAA"
# print(calculate(s+s+s+s+s+s,codon_table,AMINO_TABLE))
sequence2 = "ATAATGACAAACAAAAGTAGACTGCCGCATCAAGTGGCTGATGAAGGGTTTTATTGTTGGATAATGACAAACAAAAGTAGACTGCCGCATCAAGTGGCTGATGAAGGGTTTTATTGTTGG"
sequence = "ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"
S = "ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"
S2 = "ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"

# print(len(S2)/3)
# print(calculate(S2+"ATAATCATATGG"))
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
