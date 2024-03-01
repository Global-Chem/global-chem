#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class AndrographisPaniculata(object):

  def __init__(self):

    self.name = 'andrographis_paniculata'

  @staticmethod
  def get_smiles():

    smiles = {
      'Andrographolides': 'CC12CCC(C(C1CCC(=C)C2CC=C3C(COC3=O)O)(C)CO)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
