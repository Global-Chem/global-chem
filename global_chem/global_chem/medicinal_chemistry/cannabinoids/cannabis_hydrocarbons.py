#!/usr/bin/env python3
#
# GlobalChem - CannabisHydrocarbons
#
# ----------------------------------

class CannabisHydrocarbons(object):

  def __init__(self):

    self.name = 'cannabis_hydrocarbons'

  @staticmethod
  def get_smiles():

    smiles = {
      'n-nonane': 'CCCCCCCCC',
      'n-decane': 'CCCCCCCCCC',
      'n-undecane': 'CCCCCCCCCCC',
      'n-dodecane': 'CCCCCCCCCCCC',
      'n-tridecane': 'CCCCCCCCCCCCC',
      'd-tetradecane': 'CCCCCCCCCCCCCC',
      '3,6-dimethyl-tridecane': 'CCCCC(C)CCC(C)CC',
      'n-pentadecane': 'CCCCCCCCCCCCCCC',
      '2,6-dimethyl tetradecane': 'CCCCCCCCC(C)CCCC(C)C',
      'n-hexadecane': 'CCCCCCCCCCCCCCCC',
      'n-heptadecane': 'CCCCCCCCCCCCCCCCC',
      '2,6-dimethyl hexadecane': 'CCCCCCCCCCCC(C)CCC(C)C',
      'n-octadecane': 'CCCCCCCCCCCCCCCCCC',
      '3,6-dimethyl heptadecane': 'CCCCCCCCCCCC(C)CC(C)CC',
      '3,7-dimethyl heptadecane': 'CCCCCCCCCCC(C)CCC(C)CC',
      'n-nonadecane': 'CCCCCCCCCCCCCCCCCCC',
      '3,6-dimethyl octadecane': 'CCC(C)CCC(C)CC',
      '3,7-dimethyl octadecane': 'CC(C)CCCC(C)CC',
      'n-eicosane': 'CCCCCCCCCCCCCCCCCCCC',
      'n-heneicosane': 'CCCCCCCCCCCCCCCCCCCCC',
      '3-methyl tricosane': 'CCCCCCCCCCCCCCCCCCCCC(C)CC',
      'n-tetracosane': 'CCCCCCCCCCCCCCCCCCCCCCCC',
      '2-methyl tetracosane': 'CCCCCCCCCCCCCCCCCCCCCCC(C)C',
      'n-pentacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCC',
      'n-hexacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCCC',
      '3-methyl-pentacosane': 'CCCCCCCCCCCCCCCCCCCCCCC(C)CC',
      '2-methyl hexacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCC(C)C',
      'n-heptacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCC',
      '3-methyl heptacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCC(C)CC',
      'n-octacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      '2-methyl octacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCC(C)C',
      '9-methyl octacosane': 'CCCCCCCCCCCCCCCCCCCC(C)CCCCCCCC',
      'n-nonacosane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      '3-methyl triacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCC(C)CC',
      'n-triacotane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      '2-methyl hentriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(C)C',
      'n-hentriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      '3-methyl hentriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(C)CC',
      'n-dotriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      '2-methyl dotriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(C)C',
      'n-tritriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      'tetra-triacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      'pentatriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      'hexatriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      'heptatriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      'octatriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
      'nonatriacontane': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
