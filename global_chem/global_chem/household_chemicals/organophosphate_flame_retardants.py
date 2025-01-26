#!/usr/bin/env python3
#
# GlobalChem - Organophosphate_flame_retardants
#
# -------------------------------

class Organophosphate_flame_retardants(object):

  def __init__(self):

    self.name = 'organophosphate_flame_retardants'

  @staticmethod
  def get_smiles():

    smiles = {
'tris(1,3-dichloro-2-propyl)phosphate' : 'C(C(CCl)OP(=O)(OC(CCl)CCl)OC(CCl)CCl)Cl',
'tris(2,3-dibromopropyl) phosphate' : 'C(C(CBr)Br)OP(=O)(OCC(CBr)Br)OCC(CBr)Br',

    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {

 }

    return smarts
