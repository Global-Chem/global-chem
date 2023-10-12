#!/usr/bin/env python3
#
# GlobalChem - CannabisPigments
#
# -----------------------------

class CannabisPigments(object):

  def __init__(self):

    self.name = 'cannabis_pigments'

  @staticmethod
  def get_smiles():

    smiles = {
      'carotene': 'CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2=C(CCCC2(C)C)C)C)C',
      'zanthophylls': 'CC1=C(C(CC(C1)O)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2C(=CC(CC2(C)C)O)C)C)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
