#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class GarciniaIndica(object):

  def __init__(self):

    self.name = 'garcinia_indica'

  @staticmethod
  def get_smiles():

    smiles = {
      'hydroxycitric acid': 'C(C(=O)O)C(C(C(=O)O)O)(C(=O)O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
