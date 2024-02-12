#!/usr/bin/env python3
#
# GlobalChem - Cannabis Esters And Lactones
#
# -----------------------------------------

class CannabisEstersAndLactones(object):

  def __init__(self):

    self.name = 'cannabis_esters_and_lactones'

  @staticmethod
  def get_smiles():

    smiles = {
      'benzyl acetate': 'CC(=O)OCc1ccccc1',
      'para ethyl benzyl acetate': 'CCC1=CC=C(C=C1)COC(=O)C',
      '3-hexenyl caproate': 'CCCCCC(=O)OCCC=CCC',
      'hexyl acetate': 'CCCCCCOC(C)=O',
      'hexyl butyrate': 'CCCCCCOC(=O)CCC',
      'hexyl isobutyrate': 'CCCCCCOC(=O)C(C)C',
      'methyl acetate': 'COC(C)=O',
      'methyl linoleate': r'CCCCC/C=C\C\C=C/CCCCCCCC(=O)OC',
      'methyl palmitate': 'CCCCCCCCCCCCCCCC(=O)OC',
      'methyl salicylate': 'COC(=O)c1ccccc1O',
      'octyl caproate': 'CCCCCCCCOC(=O)CCCCC',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
