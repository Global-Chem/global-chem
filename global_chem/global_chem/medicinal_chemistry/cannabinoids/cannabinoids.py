#!/usr/bin/env python3
#
# GlobalChem - Cannabinoids
#
# -------------------------

class Cannabinoids(object):

  def __init__(self):

    self.name = 'cannabinoids'

  @staticmethod
  def get_smiles():

    smiles = {
      'cannabigerolic acid': r'CCCCCC1=CC(=C(C(=C1C(=O)O)O)C/C=C(\C)/CCC=C(C)C)O',
      'cannabigerolic acid monomethylether': r'CCCCCC1=CC(=C(C(=C1C(=O)O)O)C/C=C(\C)/CCC=C(C)C)OC',
      'cannabigerol': r'CCCCCC1=CC(=C(C(=C1)O)C/C=C(\C)/CCC=C(C)C)O',
      'cannabigerol monomethylether': r'CCCCCC1=CC(=C(C(=C1C(=O)O)O)C/C=C(\C)/CCC=C(C)C)OC',
      'cannabigerovarinic acid': 'CCCC1=CC(=C(C(=C1C(=O)O)O)CC=C(C)CCC=C(C)C)O',
      'cannabigerovarin': r'CCCC1=CC(=C(C(=C1)O)C/C=C(\C)/CCC=C(C)C)O',
      'cannabichromenic acid': 'CCCCCC1=CC2=C(C=CC(O2)(C)CCC=C(C)C)C(=C1C(=O)O)O',
      'cannabichromene': 'CCCCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
      'cannabichromevarinic acid': 'CCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
      'cannabichromevarin': 'CCCC1=CC(=C2C=CC(OC2=C1)(C)CCC=C(C)C)O',
      'cannabidiolic acid': 'CCCCCC1=CC(=C(C(=C1C(=O)O)O)C2C=C(CCC2C(=C)C)C)O',
      'cannabidiol': 'CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O',
      'cannabidiol monomethylether': 'CCCCCC1=CC(=C(C(=C1)OC)C2C=C(CCC2C(=C)C)C)O',
      'cannabidiol c4': 'CCCCC1=CC2=C(C=C1)C(OC3=CC(=CC(=C32)O)C)(C)C',
      'cannabidivarinic acid': 'CCCC1=CC(=C(C(=C1C(=O)O)O)C2C=C(CCC2C(=C)C)C)O',
      'cannabidivarin': 'CCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O',
      'cannabidiorcol': 'CC1=CC(C(CC1)C(=C)C)C2=C(C=C(C=C2O)C)O',
      'delta-9 tetrahgdrocannabinolic acid a': 'CCCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1C(=O)O)O',
      'delta-9 tetrahgdrocannabinolic acid b': 'CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1C(=O)O)(C)C)C)O',
      'delta-9 tetrahydrocannabinol': 'CCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
      'detla-9 tetrahvdrocannabinolic acid c': 'CCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1C(=O)O)O',
      'delta-9 tetrahydrocannabinol c4': 'CCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
      'delta-9 tetrahydrocannabivarinic acid': 'CCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1C(=O)O)O',
      'delta-9 tetrahydrocannabivarin': 'CCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
      'delta-9 tetrahydrocannabiorcolic acid': 'CC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1C(=O)O)O',
      'detla-9 tetrahydrocannabiorcol': 'CC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O',
      'detla-8 tetrahydrocannabinolic acid': 'CCCCCC1=CC2=C(C3CC(=CCC3C(O2)(C)C)C)C(=C1C(=O)O)O',
      'delta-8 tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3CC(=CCC3C(OC2=C1)(C)C)C)O',
      'cannabicyclolic acid': 'CCCCCC1=CC2=C(C3C4C(C3(C)C)CCC4(O2)C)C(=C1C(=O)O)O',
      'cannabicyclol': 'CCCCCC1=CC(=C2C3C4C(C3(C)C)CCC4(OC2=C1)C)O',
      'cannabicyclovarin': 'CCCC1=CC(=C2C3C4C(C3(C)C)CCC4(OC2=C1)C)O',
      'cannabielsoic acid a': 'CCCCCC1=CC2=C(C3C(CCC(C3O2)(C)O)C(=C)C)C(=C1C(=O)O)O',
      'cannabielsoic acid b': 'CCCCCC1=CC(=C2C3C(CCC(C3OC2=C1C(=O)O)(C)O)C(=C)C)O',
      'cannabielsoin': 'CCCCCC1=CC(=C2C3C(CCC(C3OC2=C1)(C)O)C(=C)C)O',
      'cannabinolic acid a': 'CCCCCC1=CC2=C(C3=C(C=CC(=C3)C)C(O2)(C)C)C(=C1C(=O)O)O',
      'cannabinol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)O',
      'cannabinol methylether': 'CCCCCC1=CC2=C(C(=C1)OC)C3=C(C=CC(=C3)C)C(O2)(C)C',
      'cannabinol c4': 'CCCCC1=CC2=C(C=C1)C(OC3=CC(=CC(=C32)O)C)(C)C',
      'cannabivarin': 'CCCC1=CC(=C2C(=C1)OC(C3=C2C=C(C=C3)C)(C)C)O',
      'cannabiorcol': 'CC1=CC2=C(C=C1)C(OC3=CC(=CC(=C32)O)C)(C)',
      'cannabinodiol': 'CCCCCC1=CC(=C(C(=C1)O)C2=C(C=CC(=C2)C)C(=C)C)O',
      'cannabinodivarin': 'CCCC1=CC(=C(C(=C1)O)C2=C(C=CC(=C2)C)C(=C)C)O',
      'cannabitriol': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C(C(CC3)(C)O)O)(C)C)O',
      '9,10-dihydroxy-delta-6a,10a-tetrahydrocannabinol': 'CCCCCC1=CC(O)=C2C(OC(C)(C)C3=C2C(O)C(C)(O)CC3)=C1',
      '10-ethoxy-9-hydroxy-delta-6a,10a-tetrahydrocannabinol': 'CCCCCC1=CC(O)=C2C(OC(C)(C)C3=C2C(OCC)C(C)(O)CC3)=C1',
      '8,9-dihydroxy-delta-tetrahydrocannabinol': 'CCCCCC1=CC(O)=C2C(OC(C)(C)C3=C2CC(C)(O)C(O)C3)=C1',
      'cannabidiolic acid tetrahydrocannabitriol ester ': 'CCCCCC1=CC(=C2C(=C1)OC(C3=C2C(C(C(C3)OC(=O)C4=C(C(=C(C=C4CCCCC)O)C5C=C(CCC5C(=C)C)C)O)C)O)(C)C)O',
      'dehydrocannabifuran': 'CCCCCC1=CC(=C2C(=C1)OC3=C(C=CC(=C23)C(=C)C)C)O',
      'cannabifuran': 'CCCCCC1=CC(=C2C(=C1)OC3=C(C=CC(=C23)C(C)C)C)O',
      'cannabichromanon': 'CCCCCC1=CC(=C2C(=C1)OC(C(C2=O)CCC(=O)C)(C)C)O',
      'cannabicitran': 'CCCCCC1=CC2=C3C4CC(CCC4C(O2)(C)C)(OC3=C1)C',
      '10-oxo-delta-6a,10a-tetrahydrocannabinol': 'CCCCCC1=CC(O)=C2C(OC(C)(C)C3=C2C(=O)C(C)CC3)=C1',
      'delta-9 6a,10a-tetrahydrocannabinol': 'CCCCCC1=CC(O)=C2C(OC(C)(C)C3=C2CC(C)(O)C(O)C3)=C1',
      '3,4,5,6-tetrahydro-7-hydroxy-alpha,alpha-2-trimethyl-9-n-propyl-2,6-methano-2H-1-benzoxocin-5-methanol': 'CCCC1=CC(=C2C3CC(CCC3C(C)(C)O)(OC2=C1)C)O',
      '6a,9,10,10a-9,10-dihydroxyhexahydrocannabinol': 'CCCCCC1=CC(=C2C3C(CCC(C3O)(C)O)C(OC2=C1)(C)C)O',
      '6a,7,10a-trihydroxy-delta-9 tetrahydrocannabinol': 'CCCCCC1=CC(=C2C3(O)C=C(CC(O)C3(O)C(OC2=C1)(C)C)C)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
