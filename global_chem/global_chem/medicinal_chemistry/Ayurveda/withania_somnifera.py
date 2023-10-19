#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class WithaniaSomnifera(object):

  def __init__(self):

    self.name = 'withania_somnifera'

  @staticmethod
  def get_smiles():

    smiles = {
      'Withanolides': 'CC1=C(C(=O)OC(C1)C(C)(C2(CCC3(C2(CCC4C3CC(C5(C4(C(=O)C=CC5)C)O)O)C)O)O)O)C',
      'Withanolide D': 'CC1=C(C(=O)OC(C1)C(C)(C2CCC3C2(CCC4C3CC5C6(C4(C(=O)C=CC6O)C)O5)C)O)C',
      'Withaferin A': 'CC1=C(C(=O)OC(C1)C(C)C2CCC3C2(CCC4C3CC5C6(C4(C(=O)C=CC6O)C)O5)C)CO',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
