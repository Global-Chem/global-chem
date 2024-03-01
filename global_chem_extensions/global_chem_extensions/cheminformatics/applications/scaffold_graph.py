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

    def __init__(self,
                 node,
                 verbose=False,
                 ):

        self.node = node
        self.verbose = verbose

    def fetch_network_data(self):

        '''

        Fetch the Network Data

        '''

        gc = GlobalChem()
        gc.build_global_chem_network()
        iupac_list = list(gc.get_node_smiles(self.node).keys())
        smiles_list = list(gc.get_node_smiles(self.node).values())

        return iupac_list, smiles_list

    def prepare_dataframe(self, iupac_list, smiles_list):

        '''

        Prepare the DataFrame

        '''

        df = pd.DataFrame()

        df['iupac_list'] = iupac_list
        df['smiles_list'] = smiles_list

        return df

    def prepare_scaffold_graph(self, df):

        '''

        Prepare the Scaffold Graph

        '''

        tree = sg.ScaffoldTree.from_dataframe(
            df,
            smiles_column='smiles_list',
            name_column='iupac_list',
            progress=True,
        )

        return tree

    def connect(self):

        '''

        Connect the Adapter

        '''

        iupac_list, smiles_list = self.fetch_network_data()
        df = self.prepare_dataframe(
            iupac_list,
            smiles_list,
        )

        graph = self.prepare_scaffold_graph(
            df
        )

        if self.verbose:
            print ("Length of SMILES List: %s" % len(smiles_list))
            print ("Length of IUPAC List: %s" % len(iupac_list))
            print ("DataFrame: %s" % df)
            print ("Scaffold Graph: %s " % graph)

        return graph