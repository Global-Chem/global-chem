#!/usr/bin/env python3
#
# GlobalChem - Ginger Terpenoids
# ------------------------------

class ThaiGingerPolysaccharides(object):

  def __init__(self):

      self.name = 'thai_ginger_polysaccharides'

  @staticmethod
  def get_smiles():

      '''
  
      Missing:
  
  
      '''
  
      smiles = {
        'fucose' : 'C[C@H]1[C@H]([C@H]([C@@H](C(O1)O)O)O)O',
        'arabinose' : 'C1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O',
        'xylose' : 'C1[C@H]([C@@H]([C@H](C(O1)O)O)O)O',
        'rhamnose' : 'C[C@H]1[C@@H]([C@H]([C@H](C(O1)O)O)O)O',
        'mannose' : 'C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O',
        'galactose' : 'C([C@@H]1[C@H]([C@@H]([C@@H](C(O1)O)O)O)O)O',
        'glucose' : 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
        'glucuronic acid' : 'C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)C(=O)O)O)O)O)O)O',
        'galacturonic acid' : '[C@@H]1([C@H]([C@H](OC([C@@H]1O)O)C(=O)O)O)O'
      }
  
  
      return smiles

  @staticmethod
  def get_smarts():

      smarts = {
  
      }
  
      return smarts
