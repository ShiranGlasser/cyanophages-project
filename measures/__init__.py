from .Nk_measure import calculate as nk_calc
from .TAI_measure import calculate as tai_calc
from .constants import MEASURES_FOLDER, FIGURES_FOLDER, CODON_TABLE
from parsers.parser import Parser
from record.helpers import get_protein_gc_number
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


""" this class saves for each gene a csv file with its measures grades
    the current measures are - Nk and TAI"""


class Measures:

    def __init__(self, records):

        if not os.path.exists(MEASURES_FOLDER):
            os.makedirs(MEASURES_FOLDER)
        if not os.path.exists(FIGURES_FOLDER):
            os.makedirs(FIGURES_FOLDER)

        all_records_nk_grades = {}
        for record in records:
            record_genes = self.get_gene_sequence(record)
            measures_df = self.build_csv(
                record.record_id, record_genes, record.seq)
            all_records_nk_grades[record.record_id] = (measures_df['Nk'][0])
            # calculates and plot the data
            self.calculate_data(
                record.record_id, record.record_data, measures_df)

        self.plot_all_records_histograms(all_records_nk_grades)

    def build_csv(self, record_id, genes_dic, seq):

        parser = Parser()
        measure = {}
        measure["all genome"] = self.get_measures(seq)
        self.all_gc = get_protein_gc_number(None, seq)
        for gene in genes_dic:

            measure[gene] = self.get_measures(
                genes_dic[gene])  # gets the measurs results

        measure_path = os.path.join(
            MEASURES_FOLDER, 'measure_' + record_id+'.csv')

        print(f"creating {record_id} measures csv file")
        # build the df to be saved in the file
        meas_attributes = parser.transpose(measure)
        # saves the results in the gene's file
        parser.to_csv(meas_attributes, measure_path)
        return meas_attributes

    def get_measures(self, seq):
        # call the calculating methods for each measure and build a dictionary of the results
        res_dic = {}
        # translated_seq = self.translate_seq(seq)[0]
        if not self.is_valid_dna(seq):
            print("error with seq {}".format(seq))
            return
        res_dic["Nk"] = nk_calc(seq.__str__())
        res_dic["TAI"] = tai_calc(seq.__str__())
        return res_dic

    def calculate_data(self, record_id, record_data, measures_df):
        normalized_Nk = self.moving_average(measures_df["Nk"][1:], 5)
        self.plot_measures_histograms(record_id, record_data, normalized_Nk)
        self.linear_regression(record_id, record_data, measures_df)

    def plot_general(self, size, x, y, x_label, y_label, title, color, save_path, to_show=None):

        plt.figure(figsize=size)
        plt.plot(x, y, 'o', color=color)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)

        plt.tight_layout()
        # Save the plot as a PNG file
        plt.savefig(save_path)
        if to_show is not None:
            # Display the plot
            plt.show()
        # Close the plot
        plt.close()

    # HELPERS :

    def get_gene_sequence(self, genome_data):
        # return a dictionary for each gene its id and sequence
        genome_df = genome_data.record_data
        rows_num = genome_df.shape[0]
        genes_dic = {}
        for i in range(0, rows_num):

            start = genome_df.iat[i, 0].__int__()
            end = genome_df.iat[i, 1].__int__()
            genes_dic[genome_df.index[i]] = genome_data.seq[start:end]
        return genes_dic

    def translate_seq(self, seq, reading_frame=None):
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

    def translate_helper(self, start_ind, seq):
        prot = ""
        stop_ind = len(seq)-(len(seq)-start_ind) % 3
        for i in range(start_ind, stop_ind, 3):
            sub_seq = seq[i] + seq[i+1] + seq[i+2]
            prot += CODON_TABLE[sub_seq]
        return prot

    def check_seq(self, seq):
        return True

    def plot_measures_histograms(self, measure_id, genome_df, measure_df):
        self.plot_nk_by_area(measure_id, genome_df, measure_df)
        self.plot_nk_by_gc(measure_id, genome_df, measure_df)

    def plot_nk_by_area(self, measure_id, genome_df, measure_df):

        folder = os.path.join(FIGURES_FOLDER, "nk by area")
        if not os.path.exists(folder):
            os.makedirs(folder)

        # Extract Nk grades and their corresponding labels
        labels = genome_df['start'].tolist()
        # Convert labels to numeric values
        labels = [int(label) for label in labels]
        # Create the plot
        self.plot_general((10, 6), labels, measure_df, 'gene position', 'Nc Value',
                          'Nc Values for Each gene', 'purple', os.path.join(folder, measure_id+'.png'))

    def plot_nk_by_gc(self, measure_id, genome_df, measure_df):

        folder = os.path.join(FIGURES_FOLDER, "nk by GC")
        if not os.path.exists(folder):
            os.makedirs(folder)
        # Extract Nk grades and their corresponding labels
        gc_values = genome_df['gc_percentage'].tolist()
        # Convert labels to numeric values
        gc_values = [float(gc_values) for gc_values in gc_values]
        plt.xlim([0, max(gc_values)])
        # Create the plot
        self.plot_general((10, 6), gc_values, measure_df, 'gc percentage', 'Nc Value',
                          'Nc Values comparing to GC percentage in a gene', 'blue', os.path.join(folder, measure_id+'.png'))

    def plot_all_records_histograms(self, all_records_nk_grades):
        
        # Extract Nk grades and their corresponding labels
        nk_values = all_records_nk_grades.values()
        labels = all_records_nk_grades.keys()

        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.scatter(labels, nk_values)
        plt.xlabel('gene id', rotation='horizontal', ha='left', va='center')
        plt.xticks(fontsize=5)
        plt.ylabel('genome Nk Value')
        plt.title('Genome Nk Values for Each record')
        plt.xticks(rotation='vertical')

        # Save the plot as a PNG file]
        fig_path = os.path.join(FIGURES_FOLDER, "all records figures")
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

        plt.savefig(os.path.join(fig_path, 'nk_for_each_record.png'))

        # Display the plot
        plt.show()
        # Close the plot
        plt.close()

    def moving_average(self, data, window_size):

        smoothed_data = []
        half_window = window_size // 2

        for i in range(len(data)):
            start = max(0, i - half_window)
            end = min(len(data), i + half_window + 1)
            window = data[start:end]
            smoothed_value = sum(window) / len(window)
            smoothed_data.append(smoothed_value)

        return smoothed_data

    def linear_regression(self, measure_id, genome_df, measure_df):

        folder = os.path.join(FIGURES_FOLDER, "linear_regression")
        if not os.path.exists(folder):
            os.makedirs(folder)
        df_fin = pd.concat(
            [measure_df[1:], genome_df['gc_percentage'], genome_df['gc_number']], axis=1)
        df_fin = df_fin.dropna()
        X = np.array(df_fin["Nk"].iloc[:-1])
        X = X.reshape(-1, 1)
        Y = np.array(df_fin["gc_percentage"].iloc[:-1])
        model = LinearRegression()

        # Fit the model to the data
        model.fit(X, Y)

        # Get the coefficient and intercept
        coef = model.coef_[0]
        intercept = model.intercept_

        # Plot the scatter plot
        plt.scatter(X, Y, color='blue', label='Actual Data')

        # Plot the linear regression line
        plt.plot(X, model.predict(X), color='red',
                 label='Linear Regression Line')

        # Add the coefficients as text on the plot
        plt.annotate(f'Coefficient: {coef:.2f}', xy=(
            X.max(), Y.min()), xytext=(X.max() - 1, Y.min() + 2), color='black')
        plt.annotate(f'Intercept: {intercept:.2f}', xy=(
            X.min(), Y.max()), xytext=(X.min() + 0.5, Y.max() - 2), color='black')

        # Set labels and legend
        plt.xlabel('Nc')
        plt.ylabel('GC')
        plt.legend()
        # Save the plot as a PNG file
        plt.savefig(os.path.join(folder, measure_id+'.png'))

        # Close the plot
        plt.close()
