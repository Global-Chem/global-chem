#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class PiperNigrum(object):

  def __init__(self):

    self.name = 'piper_nigrum'

  @staticmethod
  def get_smiles():

    smiles = {
      'Piperine': 'C1CCN(CC1)C(=O)C=CC=CC2=CC3=C(C=C2)OCO3',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
