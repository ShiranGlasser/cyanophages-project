from measures.Nk_measure import calculate as nk_calc
from measures.TAI_measure import calculate as tai_calc
from record.constants import MEASURES_FOLDER, CODON_TABLE
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
            # import pdb; pdb.set_trace()
            measure[gene] = self.get_measures(genes_dic[gene]) #get the measurs results

        measure_path = os.path.join(MEASURES_FOLDER, 'measure_'+ record_id+'.csv')

        meas_atributea = parser.transpose(measure) #build the df to be saved in the file
        parser.to_csv(meas_atributea, measure_path) #saves the results in the gene's file

        print("creating measures csv file")

    def get_measures(self, seq):
        #call the calculating methods for each measure and build a dictionary of the results
        res_dic = {}
        if not self.check_seq(seq):
            print("error with seq {}".format(seq))
            return
        res_dic["Nk"] = nk_calc(seq, CODON_TABLE)
        res_dic["TAI"] = tai_calc(seq,CODON_TABLE)
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


    def check_seq(self, seq):
        return True