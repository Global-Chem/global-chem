#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class TribulusTerrestris(object):

  def __init__(self):

    self.name = 'tribulus_terrestris'

  @staticmethod
  def get_smiles():

    smiles = {
      'Steroidal saponins': 'CC1C2C(CC3C2(CCC4C3CCC5C4(CC(C(C5)OC6C(C(C(C(O6)CO)OC7C(C(C(C(O7)C)O)O)OC8C(C(C(C(O8)C)O)O)O)O)O)O)C)C)OC19CCC(CO9)CO',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
