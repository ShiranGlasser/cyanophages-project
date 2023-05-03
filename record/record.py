import os
from collections import Counter
from typing import Type

from Bio import Entrez, SeqIO

from parsers.parser import Parser
from .helpers import get_interregions
from .constants import ENTREZ_EMAIL, GENE_BANK_FOLDER, FASTA_FOLDER, TRNA_SCAN_FOLDER, GENES, NUCLEOTIDE
# from .trnaScan import TrnaScan

fields = {'name', 'family', 'genome_size', 'gene_length_in_genome', 'percentage_gene_length_in_genome',
          'percentage_GC_in_genome', 'percentage_GC_in_genes', 'intergene_length_in_genome',
          'percentage_intergene_length_in_genome', 'percentage_GC_in_intergene'}


def assert_gb_folder():
    if not os.path.exists(GENE_BANK_FOLDER):
        os.makedirs(GENE_BANK_FOLDER)


def assert_fasta_folder():
    if not os.path.exists(FASTA_FOLDER):
        os.makedirs(FASTA_FOLDER)


def assert_sum_of_genes(dictionary):
    gene_num = dictionary["gene"]
    number_of_occurences = sum(dictionary[gene_type] for gene_type in GENES)
    assert gene_num == number_of_occurences, "gene number is: {0} and number of the genes " \
                                             "types is: {1}".format(gene_num, number_of_occurences)


def assert_percentage(number):
    assert 0 <= number <= 100


def get_protein_gc_number(seq):
    counter_g = seq.count('G')
    counter_c = seq.count('C')
    return counter_c + counter_g


def assert_number_of_intergenes_are_greater_than_genes(size_of_genes, size_of_intergene):
    assert size_of_genes <= size_of_intergene


class Record:

    def __init__(self, record_id: str, record_family: str, parser: Type[Parser]):
        print('Creating record with id: ' + record_id)
        self.record_id = record_id
        self.record_family = record_family
        self.create_genbank_file()
        self.create_fasta_file()
        # TrnaScan(record_id, record_family)
        record_content = self.get_record_content()

        self.taxonomy = record_content.annotations['taxonomy']
        self.record_data = parser.get_record_data('data\\csv\\{}.csv'.format(record_id), record_content)
        self.main_attributes = self.get_main_attributes(record_content)

        self.seq = record_content.seq.upper()

    def search(self):
        Entrez.email = ENTREZ_EMAIL
        handle = Entrez.esearch(db=NUCLEOTIDE, term=self.record_id)
        record = Entrez.read(handle)
        a = record["IdList"][0]

        record = Entrez.read(Entrez.elink(dbfrom=NUCLEOTIDE, id=a))
        print(record)

    def get_genbank_record(self):
        with Entrez.efetch(db=NUCLEOTIDE, id=self.record_id, rettype="gb", retmode="full",
                           usehistory="true", style='gbwithparts') as handle:
            list_of_records = []
            for record in SeqIO.parse(handle, "genbank"):
                list_of_records.append(record)
                print()
            return list_of_records[0]

    def get_record_content(self):
        assert os.path.exists(GENE_BANK_FOLDER + '{}.gb'.format(self.record_id))
        file_name = GENE_BANK_FOLDER + '{}.gb'.format(self.record_id)
        with open(file_name, "r") as handle:
            for i, record_gb in enumerate(SeqIO.parse(handle, "genbank")):
                return record_gb  # next(record_gb) # the last record

    def create_genbank_file(self):
        assert_gb_folder()
        # pubDateEnd = "2012/12/27"
        # pubDateStart = "2003/7/25"
        # searchTerm = f'("{pubDateStart}"[Publication Date]: "{pubDateEnd}"[Publication Date])'
        if not os.path.exists(GENE_BANK_FOLDER + '{}.gb'.format(self.record_id)):  # if the file not exists
            Entrez.email = ENTREZ_EMAIL

            with Entrez.efetch(db=NUCLEOTIDE, id=self.record_id, rettype="gbwithparts",
                               retmode="text") as handle:  # ,  term=searchTerm
                with open(GENE_BANK_FOLDER + '{}.gb'.format(self.record_id), "w") as out_handle:
                    out_handle.write(handle.read())
                print("The file: {}.gb created".format(self.record_id))

    def create_fasta_file(self):
        assert_fasta_folder()
        if not os.path.exists(FASTA_FOLDER + '{}.fasta'.format(self.record_id)):  # if the file not exists
            Entrez.email = ENTREZ_EMAIL

            with Entrez.efetch(db=NUCLEOTIDE, id=self.record_id, rettype="fasta",
                               retmode="text") as handle:  # ,  term=searchTerm
                with open(FASTA_FOLDER + '{}.fasta'.format(self.record_id), "w") as out_handle:
                    out_handle.write(handle.read().replace(".1", ""))
                print("The file: {}.fasta created".format(self.record_id))

    # def trna_scan_parser(self):
    #     # with open(TRNA_SCAN_FOLDER+'{}.fasta'.format(self.record_id), "w") as f:
    #     if (self.record_id == "FM207411.1"):
    #         with open("data\\trna_scan\\FM207411.1.fasta", "r") as f:
    #             for line in f:
    #                 line = line.rstrip()
    #                 if line.startswith('>'):
    #                     line1 = line.split()
    #                     print(line1[0])
    #                     print(line1[1])

    #                     # print(": ".join([line1[0], line1[-1]]))  # add characters you want to strip

    def get_main_attributes(self, record_content):
        genome_size = self.record_data["length"][0]

        genes_seq_len = self.record_data.loc[self.record_data['type'] == 'gene', 'length'].sum()
        gene_length_percent_in_genome = (genes_seq_len / (genome_size * 2)) * 100
        assert_percentage(gene_length_percent_in_genome)
        assert_percentage(self.record_data["gc_percentage"][0])
        GC_in_genes_number = self.record_data.loc[self.record_data['type'] == 'gene', 'gc_number'].sum()
        percentage_of_GC_in_genes = (GC_in_genes_number / genes_seq_len) * 100
        assert_percentage(percentage_of_GC_in_genes)
        intergenes, length_of_intergenes = get_interregions(record_content)
        length_of_intergenes_gc = sum(get_protein_gc_number(intergene.seq.upper()) for intergene in intergenes)
        genes_counter_dictionary = Counter(self.record_data["type"])

        main_attributes = {"name": self.record_id, "family": self.record_family, "genome_size": genome_size,
                           "gene_length_in_genome": genes_seq_len,
                           "percentage_gene_length_in_genome": gene_length_percent_in_genome,
                           "percentage_GC_in_genome": self.record_data["gc_percentage"][0],
                           "percentage_GC_in_genes": percentage_of_GC_in_genes,
                           "intergene_length_in_genome": length_of_intergenes,
                           "percentage_intergene_length_in_genome": (length_of_intergenes / (genome_size * 2)) * 100,
                           "percentage_GC_in_intergene": (length_of_intergenes_gc / length_of_intergenes) * 100}
        fields.update(genes_counter_dictionary)
        main_attributes.update(genes_counter_dictionary)

        return main_attributes
