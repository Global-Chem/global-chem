#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class  CurcumaLonga(object):

  def __init__(self):

    self.name = 'curcuma_longa'

  @staticmethod
  def get_smiles():

    smiles = {
      'Curcumin': 'COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
