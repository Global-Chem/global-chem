#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class GlycyrrhizaGlabra(object):

  def __init__(self):

    self.name = 'glycyrrhiza_glabra'

  @staticmethod
  def get_smiles():

    smiles = {
      'Glycyrrhizinic acid': 'CC1(C2CCC3(C(C2(CCC1OC4C(C(C(C(O4)C(=O)O)O)O)OC5C(C(C(C(O5)C(=O)O)O)O)O)C)C(=O)C=C6C3(CCC7(C6CC(CC7)(C)C(=O)O)C)C)C)C',
      'Lutein': 'CC1=C(C(CC(C1)O)(C)C)C=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC2C(=CC(CC2(C)C)O)C)C)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
