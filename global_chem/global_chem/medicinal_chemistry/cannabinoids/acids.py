#!/usr/bin/env python3
#
# GlobalChem - CannabisAcids
#
# --------------------------

class CannabisAcids(object):

  def __init__(self):

    self.name = 'cannabis_acids'

  @staticmethod
  def get_smiles():

    smiles = {
      'arabinonic acid': 'C(C(C(C(C(=O)O)O)O)O)O',
      'azelaic acid': 'C(CCCC(=O)O)CCCC(=O)O',
      'cinnamic acid': 'C1=CC=C(C=C1)C=CC(=O)O',
      'citric acid': 'OC(=O)CC(O)(CC(O)=O)C(O)=O',
      'glucaric acid': 'OC(C(O)C(O)C(O)=O)C(O)C(O)=O',
      'gluconic acid': 'OCC(O)C(O)C(O)C(O)C(O)=O',
      'glyceric acid': 'OCC(O)C(O)=O',
      'p-hydroxybenzoic acid': 'OC(=O)c1ccc(O)cc1',
      'p-hydroxycinnamic acid': 'OC(=O)C=Cc1ccc(O)cc1',
      'isocitric acid': 'OC(C(CC(O)=O)C(O)=O)C(O)=O',
      'malic acid': 'OC(CC(O)=O)C(O)=O',
      'malonic acid': 'OC(=O)CC(O)=O',
      '3-methyoxy-4-hydroxycinnamic acetate': 'COC1=C(C=CC(=C1)C=CC(=O)O)[O-]',
      'phosphoric acid': 'O[P](O)(O)=O',
      'pyroglutamic acid': 'OC(=O)C1CCC(=O)N1',
      'quinic acid': 'C1C(C(C(CC1(C(=O)O)O)O)O)O',
      'succinic acid': 'OC(=O)CCC(O)=O',
      'threonic acid': 'OCC(O)C(O)C(O)=O',
      'vanillic acid': 'COc1cc(ccc1O)C(O)=O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
