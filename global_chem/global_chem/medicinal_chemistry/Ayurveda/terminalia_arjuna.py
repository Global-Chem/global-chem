#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class TerminaliaArjuna(object):

  def __init__(self):

    self.name = 'terminalia_arjuna'

  @staticmethod
  def get_smiles():

    smiles = {
      'Arjunolic acid': 'CC1(CCC2(CCC3(C(=CCC4C3(CCC5C4(CC(C(C5(C)CO)O)O)C)C)C2C1)C)C(=O)O)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
