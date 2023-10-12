#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class Coleus_forskohlii(object):

  def __init__(self):

    self.name = 'Coleus_forskohlii'

  @staticmethod
  def get_smiles():

    smiles = {
      'Forskolin': 'CC(=O)OC1C(C2C(CCC(C2(C3(C1(OC(CC3=O)(C)C=C)C)O)C)O)(C)C)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
