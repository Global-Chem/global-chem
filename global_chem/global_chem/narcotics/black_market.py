#!/usr/bin/env python3
#
# GlobalChem - Black Market
#
# -------------------------

class BlackMarket(object):

  def __init__(self):

    self.name = 'black_market'

  @staticmethod
  def get_smiles():

    smiles = {
      'xylazine': 'CC1=C(C(=CC=C1)C)NC2=NCCCS2'
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
