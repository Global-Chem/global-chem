#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class SidaCordifolia(object):

  def __init__(self):

    self.name = 'sida_cordifolia'

  @staticmethod
  def get_smiles():

    smiles = {
      'Ephedrine': 'CC(C(C1=CC=CC=C1)O)NC',
      'Isoflavones': 'C1=CC(=CC=C1C2=COC3=C(C2=O)C=CC(=C3)OC4C(C(C(C(O4)CO)O)O)O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
