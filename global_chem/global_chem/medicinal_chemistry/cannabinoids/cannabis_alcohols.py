#!/usr/bin/env python3
#
# GlobalChem - CannabisAlcohols
#
# -----------------------------

class CannabisAlcohols(object):

  def __init__(self):

    self.name = 'cannabis_alcohols'

  @staticmethod
  def get_smiles():

    smiles = {
      'methanol': 'CO',
      'ethanol': 'CCO',
      'octanol-1': 'CCCCCCCC(O)',
      'octanol-3': 'CCCCCC(O)CC',
      'nonanol-1': 'CCCCCCCCC(O)',
      'hexadecanol-1': 'CCCCCCCC(O)CCCCCC',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
