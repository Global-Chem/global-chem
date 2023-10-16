#!/usr/bin/env python3
#
# GlobalChem - Ginger Cyclic Dipeptides
# -------------------------------------

class ThaiGingerCyclicDipeptides(object):

  def __init__(self):

    self.name = 'thai_ginger_cyclic_dipeptides'

  @staticmethod
  def get_smiles():

      '''
  
      Missing:
          'cyclo-(L-Val-L-Val)': '',
          'cyclo-(L-Ala-L-lle)': '',
          'cyclo-(L-Ala-L-Phe)': '',
          'cyclo-(L-Val-L-Ala)': '',
          'cyclo-(L-Leu-L-Tyr)': '',
          'cyclo-(L-Val-L-Tyr)': '',
          'cyclo-(L-Asp-OCH3-L-Phe)': '',
          'cyclo-(L-Tyr-L-lle)' : '',
          'cyclo-(L-Glu-OCH3-L-Phe)' : '',
          'cyclo-(L-Leu-L-lle)' : 'C1CCC(C1)[PbH]'
      '''
  
      smiles = {
        'cyclo-(L-Val-L-Phe}' : 'CC(C)[C@H]1C(=O)N[C@H](C(=O)N1)CC2=CC=CC=C2',
        'cyclo-(L-Val-L-Leu)' : 'CC(C)C[C@H]1C(=O)N[C@@H](C(=O)N1)C(C)C',
        'cyclo-(L-Ala-L-Leu)' : 'C[C@H]1C(=O)N[C@H](C(=O)N1)CC(C)C',
        'cyclo-(L-Phe-L-Tyr)' : 'C1=CC=C(C=C1)C[C@H]2C(=O)N[C@H](C(=O)N2)CC3=CC=C(C=C3)O',
        'cyclo-(L-Pro-L-Tyr)' : 'C1C[C@H]2C(=O)N[C@H](C(=O)N2C1)CC3=CC=C(C=C3)O',
        'cyclo-(L-Leu-L-Phe)' : 'CC(C)C[C@H]1C(=O)N[C@H](C(=O)N1)CC2=CC=CC=C2 ',
      }
  
      return smiles

  @staticmethod
  def get_smarts():

      smarts = {
  
      }
  
      return smarts
