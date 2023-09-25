#!/usr/bin/env python3
#
# GlobalChem - Ginger Phenolics
# -----------------------------

class ThaiGingerPhenolics(object):

  def __init__(self):

    self.name = 'thai_ginger_phenolics'

  @staticmethod
  def get_smiles():

      '''
  
      Missing:
  
          '4-methoxybenzyl-O-β-D-glucopyranoside': '',
          '1-O-4-carboxylphenyl-(6-O-4-hydroxybenzoyl)-β-D-glucopyranoside': '',
      '''
  
      smiles = {
        'benzoic acid' : 'C1=CC=C(C=C1)C(=O)O',
        'p-metho-xybenzoicacid' : 'COC1=CC=CC=C1C(=O)O.COC1=CC=CC=C1C(=O)O',
        'p-hydroxybenzoic acid' : 'C1=CC(=CC=C1C(=O)O)O',
        'vanillic acid' : 'COC1=C(C=CC(=C1)C(=O)O)O',
        'methyl 3,4-dihydroxybenzoate' : 'COC(=O)C1=CC(=C(C=C1)O)O',
        'methyl (2R,3S)-2,3-dihydroxy-3-(4-methoxyphenyl) propanoate' : r'COC1=CC=C(C=C1)[C@@H]([C@H](C(=O)OC)O)O',
        'ethyl (2R,3S)-2,3-dihydroxy-3-(4-methoxyphenyl) propanocate' : r'CCOC(=O)[C@@H]([C@H](C1=CC=C(C=C1)OC)O)O',
        'trans ethyl p-methoxycinnamate' : r'CCOC(=O)/C=C/C1=CC=C(C=C1)OC',
        'ferulic acid' : r'COC1=C(C=CC(=C1)/C=C/C(=O)O)O',
        'trans p-hydroxycinnamic acid' : r'C1=CC(=CC=C1/C=C/C(=O)O)O',
        'trans p-methoxycinnamic acid' : r'COC1=CC=C(C=C1)/C=C/C(=O)O',
        'trans ethyl cinnamate' : r'CCOC(=O)/C=C/C1=CC=CC=C1',
        'cis ethyl p-methoxycinnamate' : r'CCOC(=O)/C=C\C1=CC=C(C=C1)OC',
        '4-methoxy-benzyl (E)-3-(4-methoxyp-henyl) acrylate' : r'COC1=CC=C(C=C1)/C(=C\C(=O)OCC2=CC=CC=C2)/OC',
      }
  
      return smiles

  @staticmethod
  def get_smarts():
  
      smarts = {
  
      }
  
      return smarts
