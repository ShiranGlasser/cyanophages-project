import unittest
import random
from Nk_measure import calculate

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}
AMINO_TABLE = {'I': 3, 'M': 1, 'T': 4, 'N': 2, 'K': 2, 'S': 6, 'R': 6, 'L': 6, 'P': 4, 'H': 2,
               'Q': 2, 'V': 4, 'A': 4, 'D': 2, 'E': 2, 'G': 4, 'F': 2, 'Y': 2, 'END': 3, 'C': 2, 'W': 1}


class TestNkMeasure(unittest.TestCase):

    def test_single_codon_sequence(self):
        # Test case for a sequence where each amino acid is encoded by a single codon
        sequence = "ATAATGACAAACAAAAGTAGACTGCCGCATCAAGTGGCTGATGAAGGGTTTTATTGTTGGTAG"
        expected_result = 20
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)

    def test_all_codons_sequence(self):
        # Test case for a sequence where all codons are used
        sequence = "ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"
        expected_result = 61
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)

    def test_partitial_sequence(self):
        # Test case for a longer sequence
        sequence = "ATAATCATGACAAACAAAAGCAGTTCAAGAAGGCGACTACTCCTGCCACCCCACCAAGTAGTCGCAGCCGACGAAGGATTCTACTGCTGG"
        expected_result = 30
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)

    def test_partitial_sequence(self):
        # Test case for a longer sequence
        sequence = "ATAATCATGACAAACAAAAGCAGTTCAAGAAGGCGACTACTCCTGCCACCCCACCAAGTAGTCGCAGCCGACGAAGGATTCTACTGCTGG" * 10
        expected_result = 30
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)

    def shuffle_short_sequence(self):
        expected_result = 20
        sequence_lst = []        
        sequence = ""
        for amino in AMINO_TABLE:
                if amino == 'END':
                     continue
                num = random.randint(1, 10)
                codon_lst = [f for f in CODON_TABLE.keys() if CODON_TABLE[f] == amino]
                codon_ind = random.randint(0, len(codon_lst)-1)
                sequence_lst += [codon_lst[codon_ind]]*num
        random.shuffle(sequence_lst)
        for codon in sequence_lst:
            sequence+=codon
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)

    def shuffle_long_sequence(self):
        expected_result = 61
        sequence_lst = []        
        sequence = ""
        for amino in AMINO_TABLE:
                if amino == 'END':
                     continue
                num = random.randint(AMINO_TABLE[amino], 20)
                sequence_lst += [f for f in CODON_TABLE.keys() if CODON_TABLE[f] == amino]*num
        random.shuffle(sequence_lst)
        for codon in sequence_lst:
            sequence+=codon
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)
    
def test_sequence_with_stop_codons(self):
        # Test case for a sequence with stop codons
        sequence = "ATGCGCTAGCTTACCGGTCGAAATCGATCGTAGCTAGCTAGCTAGCTAGTAG"
        expected_result = 61
        result = calculate(sequence)
        self.assertEqual(result, expected_result)

def test_empty_sequence(self):
        # Test case for an empty sequence
        sequence = ""
        expected_result = 0
        result = calculate(sequence, CODON_TABLE, AMINO_TABLE)
        self.assertEqual(result, expected_result)



if __name__ == '__main__':
    unittest.main()