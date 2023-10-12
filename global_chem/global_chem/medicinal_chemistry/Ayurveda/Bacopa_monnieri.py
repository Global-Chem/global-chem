#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class Bacopa_monnieri(object):

  def __init__(self):

    self.name = 'Bacopa_monnieri'

  @staticmethod
  def get_smiles():

    smiles = {
      'Bacoside-A': 'CC(=CCCC(C)(C1C2CCC3C(C2(CC1=O)C)(CCC4C3(CCC(C4(C)C)OC5C(C(C(C(O5)CO)OC6C(C(C(CO6)O)O)O)O)O)CO)C)O)C',
      'Bacoside-B': 'CC(=CCCC(C)(C1C2CCC3C(C2(CC1=O)C)(CCC4C3(CCC(C4(C)C)OC5C(C(C(C(O5)CO)OC6C(C(C(CO6)O)O)O)O)O)CO)C)O)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
