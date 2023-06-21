#!/usr/bin/env python3
#
# GlobalChem - CannabisVitamins
#
# -----------------------------

class CannabisVitamins(object):

  def __init__(self):

    self.name = 'cannabis_vitamins'

  @staticmethod
  def get_smiles():

    smiles = {
      'vitamin k': 'CC(C)CCCC(C)CCCC(C)CCC/C(C)=C/CC1=C(C)C(=O)c2ccccc2C1=O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
