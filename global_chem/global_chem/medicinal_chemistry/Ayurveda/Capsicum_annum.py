#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class Capsicum_annum(object):

  def __init__(self):

    self.name = 'Capsicum_annum'

  @staticmethod
  def get_smiles():

    smiles = {
      'Capsaicin': 'CC(C)C=CCCCCC(=O)NCC1=CC(=C(C=C1)O)OC',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
