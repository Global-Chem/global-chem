#!/usr/bin/env python3
#
# GlobalChem - Surface_cleaners
#
# -------------------------------

class Surface_cleaners(object):

  def __init__(self):

    self.name = 'surface_cleaners'

  @staticmethod
  def get_smiles():

    smiles = {
'sodium hydroxide' : '[OH-].[Na+]',
'hydrochloric acid' : 'Cl',
'brass' : ' [Cu].[Zn].[Pb]',
'chloroxylenol' : 'CC1=CC(=CC(=C1Cl)C)O',
'1,4- dioxane' : 'C1COCCO1',

    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {

 }

    return smarts
