#!/usr/bin/env python3
#
# GlobalChem - Detergents
#
# -------------------------------

class Detergents(object):

  def __init__(self):

    self.name = 'detergents'

  @staticmethod
  def get_smiles():

    smiles = {
'diethylene glycol' : 'C(COCCO)O',
'perchloroethylene' : 'C(=C(Cl)Cl)(Cl)Cl',
'butyl cellosolve' : 'CCCCOCCO',
'phosphate edta' : ' C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O.[O-]P(=O)([O-])[O-] ',
'ammonium hydroxide' : '[NH4+].[OH-]',
'diethyl phthalate' : 'CCOC(=O)C1=CC=CC=C1C(=O)OCC',
'dimethyl phthalate' : 'COC(=O)C1=CC=CC=C1C(=O)OC',

    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {

 }

    return smarts
