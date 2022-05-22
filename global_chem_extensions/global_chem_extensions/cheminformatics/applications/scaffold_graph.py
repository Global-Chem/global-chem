#!/usr/bin/env python3
#
# GlobalChemExtensions - Scaffold Graph
#
# --------------------------------------


# Imports
# -------


import pandas as pd
import scaffoldgraph as sg

from global_chem import GlobalChem

class ScaffoldGraphAdapter(object):

    __version__ = '0.0.1'


    def __init__(self,
                 node
                 ):

        self.node = node

        self.smiles_list, self.iupac_list = self.fetch_network_data()
        self.dataframe = self.prepare_dataframe()
        self.network = self.prepare_scaffold_graph()

    def fetch_network_data(self):

        '''

        Fetch the Network Data

        '''

        gc = GlobalChem()
        gc.build_global_chem_network()
        iupac_list = list(gc.get_node_smiles(self.node).keys())
        smiles_list = list(gc.get_node_smiles(self.node).values())

        return iupac_list, smiles_list

    def prepare_dataframe(self):

        '''

        Prepare the DataFrame

        '''

        dataframe = pd.DataFrame()
        dataframe['smiles_list'] = self.smiles_list
        dataframe['iupac_list'] = self.iupac_list

        return dataframe

    def prepare_scaffold_graph(self):

        '''

        Prepare the Scaffold Graph

        '''

        network = sg.ScaffoldTree.from_dataframe(
            self.dataframe,
            smiles_column='smiles_list',
            name_column='iupac_list',
            progress=True,
        )

        return network