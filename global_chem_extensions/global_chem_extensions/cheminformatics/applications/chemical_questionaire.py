#!/usr/bin/env python3
#
# GlobalChemExtensions - Chemical Questions

# ------------------------------------------

# Imports
# -------

import pandas as pd
from rdkit import Chem
import plotly.graph_objects as go
from global_chem import GlobalChem

class ChemicalQuestionaire(object):

    '''

    Chemical Questionaire to Produce Relations

    '''

    def __init__(
      self,
      smiles_list = [],
      aromaticity=True,
      element_diversity=True,
      functional_group_diversity=True,
      stereoisomer=True,
      global_chem_node='organic_and_inorganic_bronsted_acids'
    ):

      self.smiles_list = self.test_smiles(smiles_list)
      self.aromaticity = aromaticity
      self.element_diversity = element_diversity
      self.functional_group_diversity = functional_group_diversity
      self.stereoisomer = stereoisomer
      self.global_chem_node = self.get_global_chem_patterns(global_chem_node)

    def test_smiles(self, smiles_list):

        '''

        Test the SMILES

        '''

        validate_smiles = []

        for molecule in smiles_list:

            try:
              validate_smiles.append(Chem.MolFromSmiles(molecule))
            except:
              print ('SMILES Not Accepted: %s' % molecule)

        return validate_smiles

    def get_global_chem_patterns(self, global_chem_node):

      '''

      Get the GlobalChem SMILES

      '''

      gc = GlobalChem()
      gc.build_global_chem_network()

      node_items = gc.get_node_smiles(global_chem_node)
      return node_items

    def determine_questions(self):

      '''

      Set-up the dataframe

      '''

      molecule_dataframe = pd.DataFrame()
      


    def plot_data(self):


if __name__ == '__main__':

    chemical_questionaire = ChemicalQuestionaire()



