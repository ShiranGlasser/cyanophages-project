import unittest
import random
from Nk_measure import calculate
from constants import AMINO_TABLE, CODON_TABLE

class TestNkMeasure(unittest.TestCase):

    def test_single_codon_sequence(self):
        # Test case for a sequence where each amino acid is encoded by a single codon
        sequence = "ATAATGACAAACAAAAGTAGACTGCCGCATCAAGTGGCTGATGAAGGGTTTTATTGTTGGTAG"
        expected_result = 20
        result = calculate(sequence)
        self.assertEqual(result, expected_result)

    def test_all_codons_sequence(self):
        # Test case for a sequence where all codons are used
        sequence = "ATAATCATTATGACAACCACGACTAACAATAAAAAGAGCAGTAGAAGGCTACTCCTGCTTCCACCCCCGCCTCACCATCAACAGCGACGCCGGCGTGTAGTCGTGGTTGCAGCCGCGGCTGACGATGAAGAGGGAGGCGGGGGTTCATCCTCGTCTTTCTTTTACTATTGCTGTTGGTTATTG"
        expected_result = 61
        result = calculate(sequence)
        self.assertEqual(result, expected_result)

    def test_partitial_sequence(self):
        # Test case for a longer sequence
        sequence = "ATAATCATGACAAACAAAAGCAGTTCAAGAAGGCGACTACTCCTGCCACCCCACCAAGTAGTCGCAGCCGACGAAGGATTCTACTGCTGG"
        expected_result = 30
        result = calculate(sequence)
        self.assertEqual(result, expected_result)

    def test_partitial_sequence(self):
        # Test case for a longer sequence
        sequence = "ATAATCATGACAAACAAAAGCAGTTCAAGAAGGCGACTACTCCTGCCACCCCACCAAGTAGTCGCAGCCGACGAAGGATTCTACTGCTGG" * 10
        expected_result = 30
        result = calculate(sequence)
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
        result = calculate(sequence)
        print(result)
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
        result = calculate(sequence)
        print(result)
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
        result = calculate(sequence)
        self.assertEqual(result, expected_result)



if __name__ == '__main__':
    # unittest.main()
    t=TestNkMeasure()
    # print(t.shuffle_short_sequence())
    print(t.shuffle_long_sequence())