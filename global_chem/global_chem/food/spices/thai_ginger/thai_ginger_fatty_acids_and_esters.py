#!/usr/bin/env python3
#
# GlobalChem - Ginger Fatty Acids & Esters 
# ----------------------------------------

class ThaiGingerFattyAcidsAndEsters(object):

  def __init__(self):

      self.name = 'thai_ginger_fatty_acids_and_esters'

  @staticmethod
  def get_smiles():
  
      '''
  
      Missing :
  
      '''
      smiles = {
        'stearic acid' : 'CCCCCCCCCCCCCCCCCC(=O)O',
        'dec-5-enoic acid' : 'CCCCC=CCCCC(=O)O',
        '2-tetradecenoic acid' : r'CCCCCCCCCCC/C=C/C(=O)O',
        'linolenic acid' : r'CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O',
        'linoleic acid' : r'CCCCC/C=C\C/C=C\CCCCCCCC(=O)O',
        'ethyl icosanoate' : 'CCCCCCCCCCCCCCCCCCCC(=O)OCC',
        'monopalmitin' : 'CCCCCCCCCCCCCCCC(=O)OCC(CO)O',
        'trimethyl citrate' : 'COC(=O)CC(CC(=O)OC)(C(=O)OC)O',
        '1,5-dimethyl citrate' : 'COC(=O)CC(CC(=O)OC)(C(=O)[O-])O',
        '5,6-dimethyl citrate' : 'CC1=C(C)C(=O)[O-]',
        '3-carboxyethyl-3-hydroxyglutaric acid 1,5-dimethyl ester' : 'C[As](CCC(=O)O)(O)O',
        'furan-2-carboxylc acid' : 'C1=COC(=C1)C(=O)O.C1=COC(=C1)C(=O)O',
        'dibutyl phthalate' : 'CCCCOC(=O)C1=CC=CC=C1C(=O)OCCCC',
        'pyroglutamyl-phenylalanine methyl ester' : 'C1=CC=C(C=C1)[Al](CO[As](C2=CC=CC=C2)N)N',
        'pyrogiutamyl-tyrosine methyl ester' : 'CC[As](COC)[Ge](C)(C)N',
        'phenylmethanol' : 'C1=CC=C(C=C1)CO',
      }
  
      return smiles

  @staticmethod
  def get_smarts():
  
      smarts = {
  
      }
  
      return smarts