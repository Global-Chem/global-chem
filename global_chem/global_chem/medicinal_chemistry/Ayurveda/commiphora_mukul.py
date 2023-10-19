#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class CommiphoraMukul(object):

  def __init__(self):

    self.name = 'commiphora_mukul'

  @staticmethod
  def get_smiles():

    smiles = {
      'Guggulsterone': 'CC=C1C(=O)CC2C1(CCC3C2CCC4=CC(=O)CCC34C)C',
      'E-Guggulsterone': 'CC=C1C(=O)CC2C1(CCC3C2CCC4=CC(=O)CCC34C)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
