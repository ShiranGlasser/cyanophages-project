from .Nk_measure import calculate as nk_calc
from .TAI_measure import calculate as tai_calc
from .constants import MEASURES_FOLDER, CODON_TABLE
from parsers.parser import Parser
import os
from Bio import Seq

""" this class saves for each gene a csv file with its measures grades
    the current measures are - Nk and TAI"""
    
class Measures:

    def __init__(self, records):
        
        all_records_dic = {}
        for record in records:
            record_genes = self.get_gene_sequence(record)
            self.build_csv(record.record_id, record_genes, record.seq)


    def build_csv(self, record_id, genes_dic, seq):

        if not os.path.exists(MEASURES_FOLDER):
            os.makedirs(MEASURES_FOLDER)

        parser = Parser()
        measure = {}
        measure["all genome"] = self.get_measures(seq)

        for gene in genes_dic:

            measure[gene] = self.get_measures(genes_dic[gene]) # gets the measurs results

        measure_path = os.path.join(MEASURES_FOLDER, 'measure_'+ record_id+'.csv')

        meas_atributea = parser.transpose(measure) #build the df to be saved in the file
        parser.to_csv(meas_atributea, measure_path) #saves the results in the gene's file

        print("creating measures csv file")

    def get_measures(self, seq):
        #call the calculating methods for each measure and build a dictionary of the results
        res_dic = {}

        # translated_seq = self.translate_seq(seq)[0]
        if not self.is_valid_dna(seq):
            print("error with seq {}".format(seq))
            return
        res_dic["Nk"] = nk_calc(seq.__str__())
        res_dic["TAI"] = tai_calc(seq.__str__())
        return res_dic

    def get_gene_sequence(self, genome_data):
        #return a dictionary for each gene its id and sequence
        
        genome_df = genome_data.record_data
        rows_num = genome_df.shape[0]
        genes_dic = {}
        for i in range (0, rows_num):
            
            start = genome_df.iat[i, 0].__int__()
            end = genome_df.iat[i, 1].__int__()
            genes_dic[genome_df.index[i]] = genome_data.seq[start:end]
        return genes_dic


    def translate_seq (self, seq, reading_frame=None):
        if not self.is_valid_dna(seq):
            return ""
        if reading_frame == None:
            return [self.translate_helper(0, seq), self.translate_helper(1, seq), self.translate_helper(2, seq)]
        else:
            return self.translate_helper(reading_frame-1, seq)

    def is_valid_dna(self, seq):
        str = "ACGTacgt"
        for i in range(0, len(seq)):
            if seq[i] not in str:
                return False
        return True
    
    def translate_helper (self, start_ind, seq):
        prot = ""
        stop_ind = len(seq)-(len(seq)-start_ind) % 3
        for i in range (start_ind, stop_ind, 3):
            sub_seq = seq[i] + seq[i+1] + seq[i+2]
            prot += CODON_TABLE[sub_seq]
        return prot
    
    def check_seq(self, seq):
        return True