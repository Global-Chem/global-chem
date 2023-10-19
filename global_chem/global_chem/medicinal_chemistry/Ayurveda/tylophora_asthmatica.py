#!/usr/bin/env python3
#
# GlobalChem - Ayurveda Drug Component
#
# ------------------------------

class  TylophoraAsthmatica(object):

  def __init__(self):

    self.name = 'tylophora_asthmatica'

  @staticmethod
  def get_smiles():

    smiles = {
       'Total-alkaloids': 'CN1C2CCC1CC(C2)OC(=O)C3CCC(C4=CC=CC=C34)(C5=CC=CC=C5)C(=O)OC6CC7CCC(C6)N7C',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
