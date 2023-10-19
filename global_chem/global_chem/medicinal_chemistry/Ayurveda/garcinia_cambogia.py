#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class GarciniaCambogia(object):

  def __init__(self):

    self.name = 'garcinia_cambogia'

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
