#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class PhyllanthusAmarus(object):

  def __init__(self):

    self.name = 'phyllanthus_amarus'

  @staticmethod
  def get_smiles():

    smiles = {
      'Phyllanthin': 'COCC(CC1=CC(=C(C=C1)OC)OC)C(CC2=CC(=C(C=C2)OC)OC)COC',
      'Hypophyllanthin': 'COCC1CC2=CC(=C3C(=C2C(C1COC)C4=CC(=C(C=C4)OC)OC)OCO3)OC',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
