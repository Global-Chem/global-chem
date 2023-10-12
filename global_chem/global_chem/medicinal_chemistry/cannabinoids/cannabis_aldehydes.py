#!/usr/bin/env python3
#
# GlobalChem - CannabisAldehydes
#
# ------------------------------

class CannabisAldehydes(object):

  def __init__(self):

    self.name = 'cannabis_aldehydes'

  @staticmethod
  def get_smiles():

    smiles = {
      'acetaldehyde': 'CC=O',
      'isobutyraldehyde': 'CC(C)C=O',
      'pentanal': 'CCCCC=O',
      'hexanal': 'CCCCCC=O',
      'heptanal': 'CCCCCCC=O',
      'octanal': 'CCCCCCCC=O',
      'nonanal': 'CCCCCCCCC=O',
      'decanal': 'CCCCCCCCCC=O',
      'undecanal': 'CCCCCCCCCCC=O',
      'dodecanal': 'CCCCCCCCCCCC=O',
      'tridecanal': 'CCCCCCCCCCCCC=O',
      'p-ethylbenzaldehyde': 'CCc1ccc(C=O)cc1',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
