#!/usr/bin/env python3
#
# GlobalChem - CannabisKetones
#
# ----------------------------

class CannabisKetones(object):

  def __init__(self):

    self.name = 'cannabis_ketones'

  @staticmethod
  def get_smiles():

    smiles = {
      'acetone': 'CC(C)=O',
      'heptanone-2': 'CCCCCC(=O)C',
      '2-methyl-2heptene-6-one': 'CC(C)=CCCC(C)=O',
      'decanone-2': 'CCCCCCCCC(=O)C',
      'undecanone-2': 'CCCCCCCCCC(=O)C',
      'dodecanone-2': 'CCCCCCCCCCC(=O)C',
      'pentadecanone-2': 'CCCCCCCCCCCCCC(=O)C',
      'octanone-3': 'CCCCCC(=O)CC',
      '2,2,6-trimethyl cyclohexanone': 'CC1CCCC(C)(C)C1=O',
      '2,2,6-trimethyl-5-cyclohexenone': 'CC1=CCCC(C)(C)C1=O',
      '3-decene-5-one': 'CCCCCC(=O)C=CCC',
      '6,10-dimethyl undecanone-2': 'CC(C)CCCC(C)CCCC(=O)C',
      '6,10,14-trimethyl pentadecanone-2': 'CC(C)CCCC(C)CCCC(C)CCCC(=O)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
