#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class ZingiberOfficinale(object):

  def __init__(self):

    self.name = 'zingiber_officinale'

  @staticmethod
  def get_smiles():

    smiles = {
      'Gingerols': 'CCCCCC(CC(=O)CCC1=CC(=C(C=C1)O)OC)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
