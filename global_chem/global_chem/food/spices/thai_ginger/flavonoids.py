#!/usr/bin/env python3
#
# GlobalChem - Ginger Flavonoids
# ------------------------------
class ThaiGingerFlavonoids(object):

  def __init__(self):

      self.name = 'thai_ginger_flavonoids'

  @staticmethod
  def get_smiles():

      '''
  
      Missing:
  
      '''
  
      smiles = {
        'kaempferol' : 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O',
        'luteolin': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O',
        'kaempferide' : 'COC1=CC=C(C=C1)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O'
      }
  
      return smiles

  @staticmethod
  def get_smarts():
  
      smarts = {
  
      }
  
      return smarts
