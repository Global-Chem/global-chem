#!/usr/bin/env python3
#
# GlobalChem - Cannabis Nitrogenous Compounds
#
# --------------------------------------------

class CannabisNitrogenousCompounds(object):

  def __init__(self):

    self.name = 'cannabis_nitrogenous_compounds'

  @staticmethod
  def get_smiles():

    smiles = {
      'choline': 'C[N+](C)(C)CCO',
      'trigonelline': 'C[N+]1=CC=CC(=C1)C(=O)[O-]',
      'muscarine': 'CC1C(CC(O1)C[N+](C)(C)C)O',
      'l-plus isoleucine betaine': 'C[N+](C)(C)CC(=O)[O-]',
      'neurine': 'C[N+](C)(C)C=C.[OH-]',
      'piperidine': 'C1CCNCC1',
      'hordenine': 'CN(C)CCC1=CC=C(C=C1)O',
      'ammonia': 'N',
      'methylamine': 'CN',
      'ethylamine': 'CCN',
      'n-propylamine': 'CCCN',
      'n-butylamine': 'CCCCN',
      'iso-butylamine': 'CC(C)CN',
      'secbutylamine': 'CCC(C)N',
      'dimethylamine': 'CNC',
      'pyrrolidine': 'C1CCNC1',
      'cannabisativine': 'CCCCCC(C(C1C=CCC2N1CCCNCCCCNC(=O)C2)O)O',
      'anhydrocannabisativine': 'CC1(C(C(O)C2CCCC2)O)C=CCC(C3)N1CCCNCCCCNC3=O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
