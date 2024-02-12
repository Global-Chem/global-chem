#!/usr/bin/env python3
#
# GlobalChem - Asteroid Ryugu
#
# -------------------------------

class AsteroidRyugu(object):

  def __init__(self):

    self.name = 'asteroid_ryugu'

  @staticmethod
  def get_smiles():

    smiles = {
      'methylamine': 'CN',
      'acetic_acid': 'CC(=O)O',
      'alanine': 'CC(C(=O)O)N',
      'isovaline': 'CCC(C)(C(=O)O)N',
      'naphthalene': 'C1=CC=C2C=CC=CC2=C1',
      'phenanthrene': 'C1=CC=C2C(=C1)C=CC3=CC=CC=C32',
      'pyrene': 'C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2',
      'fluoranthene': 'C1=CC=C2C(=C1)C3=CC=CC4=C3C2=CC=C4',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
