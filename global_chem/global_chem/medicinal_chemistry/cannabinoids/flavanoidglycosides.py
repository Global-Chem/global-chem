#!/usr/bin/env python3
#
# GlobalChem - CannabisFlavanoidGlycosides
#
# ----------------------------------------

class CannabisFlavanoidGlycosides(object):

  def __init__(self):

    self.name = 'cannabis_flavanoid_glycosides'

  @staticmethod
  def get_smiles():

    smiles = {
      'apigenin-7-O-para coumarylglucoside': r'OC1C(COC(=O)\C=C\C2=CC=C(O)C=C2)OC(OC2=CC(O)=C3C(=O)C=C(OC3=C2)C2=CC=C(O)C=C2)C(O)C1',
      'cosmosioside': 'OCC1OC(Oc2cc(O)c3C(=O)C=C(Oc3c2)c4ccc(O)cc4)C(O)C(O)C1O',
      'apigenin-O-glycoside': 'C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)OC4C(C(C(C(O4)CO)O)O)O)O)O',
      'isovitesin-7-O-rhamnoglucoside': 'CC1C(C(C(C(O1)OC2=C(C(=C3C(=C2)OC(=CC3=O)C4=CC=C(C=C4)O)O)C5C(C(C(C(O5)CO)O)O)O)O)O)O',
      'kaempferol-O-glycoside': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)OC4C(C(C(C(O4)CO)O)O)O)O)O)O',
      'luteolin-O-glycoside': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)OC4C(C(C(C(O4)CO)O)O)O',
      'orientin': 'OCC1OC(C(O)C(O)C1O)c2c(O)cc(O)c3C(=O)C=C(Oc23)c4ccc(O)c(O)c4',
      'orientin-O-glucoside': 'C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)C4C(C(C(C(O4)CO)O)O)O)O)O',
      'orientin-7-O-rhamnoglucoside': 'CC1C(C(C(C(O1)OC2=C(C3=C(C(=C2)O)C(=O)C=C(O3)C4=CC(=C(C=C4)O)O)C5C(C(C(C(O5)CO)O)O)O)O)O)O',
      'quercetin-O-glucoside': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)OC4C(C(C(C(O4)CO)O)O)O)O)O)O)O',
      'vitexin-7-O-g"-glucoside': 'C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)OC4C(C(C(C(O4)CO)O)O)O)C5C(C(C(C(O5)CO)O)O)O)O',
      'vitexin-O-glucoside': 'C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)C4C(C(C(C(O4)CO)O)O)O)OC5C(C(C(C(O5)CO)O)O)O',
      'vitexin-O-rhamnoglucoside': 'CC1C(C(C(C(O1)OC2C(C(C(OC2C3=C(C=C(C4=C3OC(=CC4=O)C5=CC=C(C=C5)O)O)O)CO)O)O)O)O)O',
      '2"-O-glucopyranosylvitexin': 'C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)C4C(C(C(C(O4)CO)O)O)OC5C(C(C(C(O5)CO)O)O)O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
