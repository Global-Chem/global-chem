#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class PicrorhizaKurroa(object):

  def __init__(self):

    self.name = 'picrorhiza_kurroa'

  @staticmethod
  def get_smiles():

    smiles = {
      'Kutkin': 'COC1=C(C=CC(=C1)C(=O)OC2C(C(C(C(O2)CO)O)O)O)OC(=O)C=CC3=CC=CC=C3.O.O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
