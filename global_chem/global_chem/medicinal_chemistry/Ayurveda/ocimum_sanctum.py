#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class  OcimumSanctum(object):

  def __init__(self):

    self.name = 'ocimum_sanctum'

  @staticmethod
  def get_smiles():

    smiles = {
      'Ursolic acid': 'CC1CCC2(CCC3(C(=CCC4C3(CCC5C4(CCC(C5(C)C)O)C)C)C2C1C)C)C(=O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
