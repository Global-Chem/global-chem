#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class MucunaPruriens(object):

  def __init__(self):

    self.name = 'mucuna_pruriens'

  @staticmethod
  def get_smiles():

    smiles = {
      'L-Dopa': 'C1=CC(=C(C=C1CC(C(=O)O)N)O)O',
      'Catecholamines': 'C1=CC(=C(C(=C1)O)O)N',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
