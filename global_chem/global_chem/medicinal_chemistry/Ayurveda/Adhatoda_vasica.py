#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class Adhatoda_vasica(object):

  def __init__(self):

    self.name = 'Adhatoda_vasica'

  @staticmethod
  def get_smiles():

    smiles = {
      'Vasicine': 'C1CN2CC3=CC=CC=C3N=C2[C@H]1O',
      'Total-alkaloids': 'CN1C2CCC1CC(C2)OC(=O)C3CCC(C4=CC=CC=C34)(C5=CC=CC=C5)C(=O)OC6CC7CCC(C6)N7C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
