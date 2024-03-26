#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class CamelliaSinensis(object):

  def __init__(self):

    self.name = 'camellia_sinensis'

  @staticmethod
  def get_smiles():

    smiles = {
      'Catechins': 'C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O',
      'Caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
