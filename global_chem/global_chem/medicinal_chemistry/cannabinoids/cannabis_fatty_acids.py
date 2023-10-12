#!/usr/bin/env python3
#
# GlobalChem - Cannabis Fatty Acids
#
# ---------------------------------

class CannabisFattyAcids(object):

  def __init__(self):

    self.name = 'cannabis_fatty_acids'

  @staticmethod
  def get_smiles():

    smiles = {
      'arachidic acid': 'CCCCCCCCCCCCCCCCCCCC(O)=O',
      'behenic acid': 'CCCCCCCCCCCCCCCCCCCCCC(O)=O',
      'eicosadienic acid': 'CCCCCCCCCCCCCCCC=CC=CC(O)=O',
      'eicosemic acid': r'O=C(O)CCCCCCC\C=C/CCCCCCCCCC',
      'linoleic': r'CCCCC/C=C/C/C=C/CCCCCCCC(O)=O',
      'linolenic acid': r'CC/C=C/C/C=C/C/C=C/CCCCCCCC(O)=O',
      'myristic acid': 'CCCCCCCCCCCCCC(O)=O',
      'oleic acid': r'CCCCCCCC\C=C/CCCCCCCC(O)=O',
      'palmitic acid': 'CCCCCCCCCCCCCCCC(O)=O',
      'palmitoleic acid': r'CCCCCC\C=C/CCCCCCCC(O)=O',
      'sativic acid': 'CCCCCC(O)C(O)CC(O)C(O)CCCCCCCC(O)=O',
      'stearic acid': 'CCCCCCCCCCCCCCCCCC(O)=O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
