#!/usr/bin/env python3
#
# GlobalChem - Ginger Diaryl Heptanoids
# -------------------------------------

class ThaiGingerDiarylHeptanoids(object):

  def __init__(self):

    self.name = 'thai_ginger_diaryl_heptanoids'

  @staticmethod
  def get_smiles():
  
      '''
  
      Missing:
          'hedycoropyran B' : '',
      '''
  
      smiles = {
        '(3R,5S)-3,5-dihydroxy-1,7-bis(3,4- dihydroxyphenyl) heptane' : 'C1=CC(=C(C=C1CC[C@H](C[C@H](CCC2=CC(=C(C=C2)O)O)O)O)O)O',
        '(1R,3R,5R)-1,5-epoxy-3-hydroxy-1-(3,4-dihydroxyphenyl)-7-(3,4-dihydroxypheny) heptane' : 'C(C[C@H](C[C@H](C[C@H](C1=CC=C(C(=C1)O)O)OO)O)O)C2=CC=C(C(=C2)O)O',
        '(1R,3R,5R)-1,5-epoxy-3-hydroxy-1-(3,4-dihydroxyphenyl)-7-(4-hydroxyphenyl) heptane 3-O-p-D-glucopyranoside' : 'C1=C(C=CC(=C1)O)CC[C@H](C[C@H](C[C@H](C2=CC=C(C(=C2)O)O)O[Hg])O)O[67O]',
        'phaeoheptanoxide' : 'CCCCCCC(=O)[O-]',
        '1-(4-hydroxy-3-methoxyphenyl)-7-(4-hydroxyphenyl) heptane-1,2,3,5,6-pentanol' : 'COC1=CC(=CC=C1O)C(C(C(CC(C(CC2=CC=C(C=C2)O)O)O)O)O)O',
        'kaempsulfonic acid A': 'OS(=O)(=O)[SH+]',
      }
  
      return smiles

  @staticmethod
  def get_smarts():

      smarts = {
  
      }
  
      return smarts