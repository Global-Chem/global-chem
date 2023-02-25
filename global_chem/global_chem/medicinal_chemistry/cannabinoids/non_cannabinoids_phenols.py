#!/usr/bin/env python3
#
# GlobalChem - Non-Cannabinoid Phenols
#
# -------------------------------------

class NonCannabinoidPhenols(object):

  def __init__(self):

    self.name = 'non_cannabinoid_phenols'

  @staticmethod
  def get_smiles():

    smiles = {
      'acetylcannabispirol': 'COc1cc2CCC3(CCC(CC3)OC(=O)C)c2c(O)c1',
      'cannabispiradienone': 'COC1=CC2=C(C(=C1)O)C3(CC2)C=CC(=O)C=C3',
      'beta cannabispiranol': 'COC1=CC2=C(C(=C1)O)C3(CCC(CC3)O)CC2',
      'cannabispirenone': 'COC1=CC2=C(C(=C1)O)C3(CCC(=O)C=C3)CC2',
      'cannabispirenone-isomer': 'COC1=CC2=C(C(=C1)O)C3(CCC(=O)C=C3)CC2',
      'cannabispirone': 'COC1=CC2=C(C(=C1)O)C3(CCC(=O)CC3)CC2',
      '3-[2-(4-hydroxyphenyl)ethyl]-5-methoxyphenol': 'COc1cc(O)cc(CCc2ccc(O)cc2)c1',
      '3-[2-(3-hydroxy-4-methoxyphenyl)ethyl]-5-methoxyphenol': 'COc1cc(O)cc(CCc2ccc(OC)c(O)c2)c1',
      '3-[2-(3-isoprenyl-4-hydroxy-5-methoxy-phenyl)ethyl]-5-methoxyphenol': 'COc1cc(O)cc(CCc2cc(OC)c(O)c(C=CC(C)=C)c2)c1',
      'canniprene': 'CC(=CCC1=C(C=CC(=C1O)OC)CCC2=CC(=CC(=C2)OC)O)C',
      'eugenol': 'COc1cc(CC=C)ccc1O',
      'isoeugenol': 'COc1cc(C=CC)ccc1O',
      'anethol': 'COc1ccc(C=CC)cc1',
      'methyleugenol': 'COc1ccc(CC=C)cc1OC',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
