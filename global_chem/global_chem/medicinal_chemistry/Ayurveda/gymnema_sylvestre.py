#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class GymnemaSylvestre(object):

  def __init__(self):

    self.name = 'gymnema_sylvestre'

  @staticmethod
  def get_smiles():

    smiles = {
      'Gymnemic acid': 'CCC(C)C(=O)OC1C(C2(C(CC1(C)C)C3=CCC4C5(CCC(C(C5CCC4(C3(CC2O)C)C)(C)CO)OC6C(C(C(C(O6)C(=O)O)O)OC7C(=O)C(C(C(O7)CO)O)O)O)C)CO)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
