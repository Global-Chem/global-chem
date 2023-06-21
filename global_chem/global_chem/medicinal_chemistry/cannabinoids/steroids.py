#!/usr/bin/env python3
#
# GlobalChem - CannabisSteroids
#
# --------------------------------

class CannabisSteroids(object):

  def __init__(self):

    self.name = 'cannabis_steroids'

  @staticmethod
  def get_smiles():

    smiles = {
      'campesterol': 'CC(C)C(C)CCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C',
      'campest-5-en-3beta-ol-7-one': 'CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2=CCC4C3(CCC(C4)O)C)C',
      'ergosterol': r'CC(C)C(C)/C=C/C(C)C1CCC2C3=CC=C4CC(O)CCC4(C)C3CCC12C',
      'beta-sitosterol': 'CCC(CCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C)C(C)C',
      '5alpha-stigmasta-7,24-dien-3beta-ol': 'CCC(CCC(C)C1CCC2C3=CCC4CC(O)CCC4(C)C3CCC12C)=C(C)C',
      'stigmasta-5,22-dien-3beta-ol-7-one': 'CCC(C=CC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)OC(=O)C)C)C)C(C)C',
      'stigmast-5-en-3beta-ol-7-one': 'CCC(CCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C)C(C)C',
      'stigmast-4-en-3-one': 'CCC(CCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C)C(C)C',
      'stigmasterol': r'CCC(/C=C/C(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C)C(C)C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
