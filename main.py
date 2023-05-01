from typing import Tuple

from parsers.helpers import pd
from figures import figures
from figures.figures import Figures
from record.record import Record
from records.records import Records
from typing import Type
from measures import Measures as ms


def species_old() -> dict[str, str]:
    return {
        # "KM034562.1":"KM034562.1",
        # "NC_045512.2": "NC_045512.2", # cororna virus
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "P_SSM4": "NC_006884.1",  # circular PHG 26-MAR-2010
        "P_SSM2": "NC_006883.1",  # circular PHG 26-MAR-2010
        "Syn9": "NC_008296.1",  # circular PHG 28-NOV-2007
        "S_RSM4": "FM207411.1",  # circular PHG 22-SEP-2009
        "S_PM2": "NC_006820.1",  # circular PHG 11-OCT-2021 !!!!!!!!!!!!!
        # ------Podoviruses:------
        "P_SSP7": "NC_006882.1",  # circular PHG 19-NOV-2010
        "Syn5": "NC_009531.1",  # linear   PHG 20-DEC-2020 !!!!!!!!!!!!!
        "P60": "NC_003390.1",  # linear   PHG 17-APR-2009

        # '''***Cyanobacteria:***'''
        # ------Prochlorococcus:------
        # - - - - - - HL-Prochlorococcus:
        "MED4": "NC_005072.1",  # circular CON 10-OCT-2021  !!!!!!!!!!!!!
        "MIT_9515": "NC_008817.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        # "MIT_9312": "NC_007577.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9215": "NC_009840.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        # - - - - - - LL-Prochlorococcus:
        "NATL1A": "NC_008819.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "NATL2A": "NC_007335.1",  # circular BCT 29-NOV-2007
        "SS120": "NC_005042.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        # "MIT_9211": "NC_009976.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "MIT_9303": "NC_008820.1",  # circular CON 25-NOV-2016 !!!!!!!!!!!!!
        # "MIT_9313": "NC_005071.1",  # circular CON 13-DEC-2020 !!!!!!!!!!!!!

        # ------Synechococcus:------
        "WH_8102": "NC_005070.1",  # circular CON 10-OCT-2021 !!!!!!!!!!!!!
        "WH_7803": "NC_009481.1",  # circular CON 17-APR-2017 !!!!!!!!!!!!!

    }


def species() -> dict[str, str]:
    print('Creting species.')
    return {
        # "KM034562.1":"KM034562.1",
        # "NC_045512.2": "NC_045512.2", # cororna virus
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "NC_006884.1": "Myoviridae",  # [Phage: Myoviridae P-SSM4], circular PHG 26-MAR-2010
        "NC_006883.1": "Myoviridae",  # [Phage: Myoviridae P_SSM2], circular PHG 26-MAR-2010
        "NC_008296.1": "Myoviridae",  # [Phage: Myoviridae Syn9], circular PHG 28-NOV-2007
        "FM207411.1": "Myoviridae",  # [Phage: Myoviridae S_RSM4], circular PHG 22-SEP-2009
        "NC_006820.1": "Myoviridae",  # [Phage: Myoviridae S_PM2], circular PHG 11-OCT-2021 !
        # # # ------Podoviruses:------
        "NC_006882.1": "Podoviridae",  # [Phage: Podoviridae P_SSP7], circular PHG 19-NOV-2010
        "NC_009531.1": "Podoviridae",  # [Phage: Podoviridae Syn5], linear   PHG 20-DEC-2020 !
        "NC_003390.1": "Podoviridae",  # [Phage: Podoviridae P60], linear   PHG 17-APR-2009

        # '''***Cyanobacteria:***'''
        # ------Prochlorococcus:------
        # - - - - - - HL-Prochlorococcus:
        "NC_005072.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MED4], circular CON 10-OCT-2021  !
        "NC_008817.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9515], circular CON 10-OCT-2021 !
        "NC_007577.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9312], circular CON 10-OCT-2021 !
        "NC_009840.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9215], circular CON 10-OCT-2021 !
        # # - - - - - - LL-Prochlorococcus:
        "NC_008819.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus NATL1A], circular CON 10-OCT-2021 !
        "NC_007335.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus NATL2A], circular BCT 29-NOV-2007
        "NC_005042.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus SS120], circular CON 10-OCT-2021 !
        "NC_009976.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9211],
        # circular CON 10-OCT-2021 !
        "NC_008820.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9303], circular CON 25-NOV-2016 !
        "NC_005071.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9313],
        # circular CON 13-DEC-2020 !

        # ------Synechococcus:------
        "NC_005070.1": "Synechococcus",  # [Cyanobacteria: Synechococcus WH_8102], circular CON 10-OCT-2021 !
        "NC_009481.1": "Synechococcus",  # [Cyanobacteria: Synechococcus WH_7803], circular CON 17-APR-2017 !
    }


def species_updated() -> dict[str, str]:
    #print('Creting species.')
    return {
        # "KM034562.1":"KM034562.1",
        # "NC_045512.2": "NC_045512.2", # cororna virus
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "NC_006884.1": "Myoviridae",  # [Phage: Myoviridae P-SSM4], circular PHG 26-MAR-2010
        "NC_006883.1": "Myoviridae",  # [Phage: Myoviridae P_SSM2], circular PHG 26-MAR-2010
        "NC_008296.1": "Myoviridae",  # [Phage: Myoviridae Syn9], circular PHG 28-NOV-2007
        "FM207411.1": "Myoviridae",  # [Phage: Myoviridae S_RSM4], circular PHG 22-SEP-2009
        "NC_006820.1": "Myoviridae",  # [Phage: Myoviridae S_PM2], circular PHG 11-OCT-2021 !

        "NC_015280.1": "Myoviridae",  # [Phage: Myoviridae P-HM1]
        "NC_015284.1": "Myoviridae",  # [Phage: Myoviridae P-HM2]
        "NC_015283.1": "Myoviridae",  # [Phage: Myoviridae P-RSM4]
        "NC_015285.1": "Myoviridae",  # [Phage: Myoviridae Syn33]#############
        "NC_015290.1": "Myoviridae",  # [Phage: Myoviridae P-SSM7]
        "NC_015282.1": "Myoviridae",  # [Phage: Myoviridae S-SM1]
        "NC_015279.1": "Myoviridae",  # [Phage: Myoviridae S-SM2]
        "NC_015289.1": "Myoviridae",  # [Phage: Myoviridae S-SSM5]
        "NC_015281.1": "Myoviridae",  # [Phage: Myoviridae S-ShM2]
        "NC_015286.1": "Myoviridae",  # [Phage: Myoviridae Syn19]

        "NC_015288.1": "Myoviridae",  # [Phage: Myoviridae Syn1]#new!
        "NC_015287.1": "Myoviridae",  # [Phage: Myoviridae S-SSM7]#new!

        # # # ------Podoviruses:------
        "NC_006882.1": "Podoviridae",  # [Phage: Podoviridae P_SSP7], circular PHG 19-NOV-2010
        "NC_009531.1": "Podoviridae",  # [Phage: Podoviridae Syn5], linear   PHG 20-DEC-2020 !
        "NC_003390.1": "Podoviridae",  # [Phage: Podoviridae P60], linear   PHG 17-APR-2009

        # '''***Cyanobacteria:***'''
        # ------Prochlorococcus:------
        # - - - - - - HL-Prochlorococcus:
        "NC_005072.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MED4], circular CON 10-OCT-2021  !
        "NC_008817.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9515], circular CON 10-OCT-2021 !
        "NC_007577.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9312], circular CON 10-OCT-2021 !
        "NC_009840.1": "HL-Prochlorococcus",  # [Cyanobacteria: HL-Prochlorococcus MIT_9215], circular CON 10-OCT-2021 !
        # # - - - - - - LL-Prochlorococcus:
        "NC_008819.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus NATL1A], circular CON 10-OCT-2021 !
        "NC_007335.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus NATL2A], circular BCT 29-NOV-2007
        "NC_005042.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus SS120], circular CON 10-OCT-2021 !
        "NC_009976.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9211],
        # circular CON 10-OCT-2021 !
        "NC_008820.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9303], circular CON 25-NOV-2016 !
        "NC_005071.1": "LL-Prochlorococcus",  # [Cyanobacteria: LL-Prochlorococcus MIT_9313],
        # circular CON 13-DEC-2020 !

        # ------Synechococcus:------
        "NC_005070.1": "Synechococcus",  # [Cyanobacteria: Synechococcus WH_8102], circular CON 10-OCT-2021 !
        "NC_009481.1": "Synechococcus",  # [Cyanobacteria: Synechococcus WH_7803], circular CON 17-APR-2017 !
    }


def viruses_and_hosts() -> dict[str, list[str]]:
    print('Creting viruses and hosts.')
    return {
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "NC_006884.1": ['LL-Prochlorococcus', 'HL-Prochlorococcus'],
        # [Phage: Myoviridae P-SSM4], circular PHG 26-MAR-2010
        "NC_006883.1": ['LL-Prochlorococcus'],  # [Phage: Myoviridae P_SSM2], circular PHG 26-MAR-2010
        "NC_008296.1": ['Synechococcus', 'LL-Prochlorococcus', 'HL-Prochlorococcus'],
        # [Phage: Myoviridae Syn9], circular PHG 28-NOV-2007
        "FM207411.1": ['Synechococcus'],  # [Phage: Myoviridae S_RSM4], circular PHG 22-SEP-2009
        "NC_006820.1": ['Synechococcus'],  # [Phage: Myoviridae S_PM2], circular PHG 11-OCT-2021 !
        # ------Podoviruses:------
        "NC_006882.1": ['HL-Prochlorococcus'],  # [Phage: Podoviridae P_SSP7], circular PHG 19-NOV-2010
        "NC_009531.1": ['Synechococcus'],  # [Phage: Podoviridae Syn5], linear   PHG 20-DEC-2020 !
        "NC_003390.1": ['Synechococcus'],  # [Phage: Podoviridae P60], linear   PHG 17-APR-2009
    }


def viruses_and_hosts_updated() -> dict[str, list[str]]:

    return {
        # '''***Cyanophages:***'''
        # ------Moyviruses:------
        "NC_006884.1": ['LL-Prochlorococcus', 'HL-Prochlorococcus'],
        # [Phage: Myoviridae P-SSM4], circular PHG 26-MAR-2010
        "NC_006883.1": ['LL-Prochlorococcus'],  # [Phage: Myoviridae P_SSM2], circular PHG 26-MAR-2010
        "NC_008296.1": ['Synechococcus', 'LL-Prochlorococcus', 'HL-Prochlorococcus'],
        # [Phage: Myoviridae Syn9], circular PHG 28-NOV-2007
        "FM207411.1": ['Synechococcus'],  # [Phage: Myoviridae S_RSM4], circular PHG 22-SEP-2009
        "NC_006820.1": ['Synechococcus'],  # [Phage: Myoviridae S_PM2], circular PHG 11-OCT-2021 !

        "NC_015280.1": ['Prochlorococcus'],  # [Phage: Myoviridae P-HM1]
        "NC_015284.1": ['Prochlorococcus'],  # [Phage: Myoviridae P-HM2]
        "NC_015283.1": ['Prochlorococcus'],  # [Phage: Myoviridae P-RSM4]
        "NC_015285.1": ['Synechococcus', 'Prochlorococcus'],  # [Phage: Myoviridae Syn33]#################
        "NC_015290.1": ['Prochlorococcus'],  # [Phage: Myoviridae P-SSM7]
        "NC_015282.1": ['Synechococcus'],  # [Phage: Myoviridae S-SM1]
        "NC_015279.1": ['Synechococcus'],  # [Phage: Myoviridae S-SM2]
        "NC_015289.1": ['Synechococcus'],  # [Phage: Myoviridae S-SSM5]
        "NC_015281.1": ['Synechococcus'],  # [Phage: Myoviridae S-ShM2]
        "NC_015286.1": ['Synechococcus'],  # [Phage: Myoviridae Syn19]

        "NC_015288.1": ['Prochlorococcus'],  # [Phage: Myoviridae Syn1]#new!
        "NC_015287.1": ['Synechococcus'],  # [Phage: Myoviridae S-SSM7]#new!

        # ------Podoviruses:------
        "NC_006882.1": ['HL-Prochlorococcus'],  # [Phage: Podoviridae P_SSP7], circular PHG 19-NOV-2010
        "NC_009531.1": ['Synechococcus'],  # [Phage: Podoviridae Syn5], linear   PHG 20-DEC-2020 !
        "NC_003390.1": ['Synechococcus'],  # [Phage: Podoviridae P60], linear   PHG 17-APR-2009
    }


def get_records(species: dict[str, str]) -> Tuple[list[Record], pd.DataFrame]:
    print('Create records of species.')
    r = Records(species)
    return r.records, r.attributes


def manage_figures(types:list[str], records:list[Record], attributes:pd.DataFrame, viruses_and_hosts:dict[str,list[str]]):
    frequencies = {}
    #TO_DEL#types = ["Podoviridae", "Myoviridae", "Prochlorococcus", "Synechococcus"]
    figures = Figures()
    x = figures.get_codon_table(11)
    for record in records:
        print(record.record_id)
        # ms(record.seq, record.record_id)
        # print("^^^^^^^^^",record.seq)
        # import pdb; pdb.set_trace()

        # start_codon = figures.get_start_codon(record.seq, x)
        # record contains list of all the genes
        # #1. get list of the CDS sequences and there start codons
        #  2. get frequency of codons  a) For cds in CDS{
        # add to frequency dict the cds frequency of the seq starting from the start codon
        #                                #                                   }
        cds_list = record.record_data.loc[record.record_data['type'] == 'CDS', "seq"].values.tolist()
        start_codon_list = record.record_data.loc[record.record_data['type'] == 'CDS', "start_codon"]
        codons = figures.get_frequency_of_codons(cds_list, start_codon_list)
        frequencies[record.record_id] = codons

        #figures.show_relative_position(record.record_data, attributes, record.record_family, record.record_id)

    mean_and_std_of_types = figures.get_mean_and_std(records, types)

    #figures.scatter_plot(attributes, viruses_and_hosts)
    #figures.bar_chart_histogram(mean_and_std_of_types, 'Genome size',
    #                             '(B) Average genome sizes of Podoviruses,'
    #                            ' Myoviruses, Prochlorococcus and Synechococcus',
    #                             'figure1_B_bar_plot_with_error_bars.png')  # figure1, B
    figures.test_stripchart(
        ["Myoviridae+Myoviridae", "Podoviridae+Podoviridae", "HL-Prochlorococcus+HL-Prochlorococcus",
          "LL-Prochlorococcus+LL-Prochlorococcus", "Synechococcus+Synechococcus"], frequencies, attributes,
        viruses_and_hosts, x.stop_codons)


if __name__ == '__main__':
    species_id_and_type = species_updated()
    # viruses_and_hosts = viruses_and_hosts_updated()
    records, attributes = get_records(species_id_and_type)
    ms( records)
    # import pdb; pdb.set_trace()

    # types = ["Podoviridae", "Myoviridae", "Prochlorococcus", "Synechococcus"]

    # manage_figures(types, records, attributes, viruses_and_hosts)
