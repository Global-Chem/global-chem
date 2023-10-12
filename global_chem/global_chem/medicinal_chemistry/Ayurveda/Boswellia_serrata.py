#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class Boswellia_serrata(object):

  def __init__(self):

    self.name = 'Boswellia_serrata'

  @staticmethod
  def get_smiles():

    smiles = {
      'Organic_acids': 'CC(C(=O)O)O',
      'Boswellic_acid': 'CC1CCC2(CCC3(C(=CCC4C3(CCC5C4(CCC(C5(C)C(=O)O)O)C)C)C2C1C)C)C',
      'Sennoside': 'C1=CC2=C(C(=C1)OC3C(C(C(C(O3)CO)O)O)O)C(=O)C4=C(C2C5C6=C(C(=CC=C6)OC7C(C(C(C(O7)CO)O)O)O)C(=O)C8=C5C=C(C=C8O)C(=O)O)C=C(C=C4O)C(=O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
