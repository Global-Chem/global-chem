#!/usr/bin/env python3
#
# GlobalChem - Bleach
#
# -------------------------------

class Bleach(object):

  def __init__(self):

    self.name = 'bleach'

  @staticmethod
  def get_smiles():

    smiles = {
'chlorine' : 'ClCl',
'sodium hypochlorite' : '[O-]Cl.[Na+]',
'hydrogen peroxide' : 'OO',
'potassium permanganate' : '[O-]Mn(=O)=O.[K+]',

    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {

 }

    return smarts
