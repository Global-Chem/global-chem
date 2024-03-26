#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class CentellaAsiatica(object):

  def __init__(self):

    self.name = 'centella_asiatica'

  @staticmethod
  def get_smiles():

    smiles = {
      'triterpenes-saponin': 'CC1C(C(C(C(O1)OC2C(C(C(OC2C(=O)O)OC3CCC4(C(C3(C)C=O)CCC5(C4CC=C6C5(CCC7(C6CC(CC7)(C)C)C(=O)OC8C(C(C(C(O8)C)OC9C(C(C(CO9)O)O)O)OC1C(C(C(C(O1)C)OC1C(C(C(C(O1)CO)O)O)O)O)OC1C(C(C(C(O1)CO)O)O)O)O)C)C)C)O)OC1C(C(C(CO1)OC1C(C(C(CO1)O)O)O)O)OC1C(C(C(C(O1)CO)O)O)O)O)O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
