#!/usr/bin/env python3
#
# GlobalChem - Cosmetics
#
# -------------------------------

class Cosmetics(object):

  def __init__(self):

    self.name = 'cosmetics'

  @staticmethod
  def get_smiles():

    smiles = {
'oxybenzone' : 'COC1=CC(=C(C=C1)C(=O)C2=CC=CC=C2)O',
'hydroquinone' : 'C1=CC(=CC=C1O)O',
'diethyl phthalate' : 'CCOC(=O)C1=CC=CC=C1C(=O)OCC',

    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {

 }

    return smarts
