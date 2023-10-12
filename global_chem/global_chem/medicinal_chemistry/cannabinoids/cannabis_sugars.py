#!/usr/bin/env python3
#
# GlobalChem - Cannabis Sugar
#
# ---------------------------

class CannabisSugars(object):

  def __init__(self):

    self.name = 'cannabis_sugars'

  @staticmethod
  def get_smiles():

    smiles = {
      'arabinose': 'C1C(C(C(C(O1)O)O)O)O',
      'fructose': 'C1C(C(C(C(O1)(CO)O)O)O)O',
      'galactose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
      'galacturonic acid': 'C1(C(C(OC(C1O)O)C(=O)O)O)O',
      'alpha and beta d-glucose': 'C(C(C(C(C(C=O)O)O)O)O)O',
      'altro-heptulose': 'C(C1C(C(C(C(O1)(CO)O)O)O)O)O',
      'd-manno-heptulose': 'C(C(C(C(C(C(=O)CO)O)O)O)O)O',
      'mannose': 'C(C1C(C(C(C(O1)O)O)O)O)O',
      'd-glycerol-d-manno-octulose': 'C(C(C(C(C(C(C(=O)CO)O)O)O)O)O)O',
      'rhamnose': 'CC1C(C(C(C(O1)O)O)O)O',
      'ribose': 'C1C(C(C(C(O1)O)O)O)O',
      'xylose': 'C1C(C(C(C(O1)O)O)O)O',
      'sucrose': 'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O',
      'maltose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
      'raffinose': 'C(C1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3(C(C(C(O3)CO)O)O)CO)O)O)O)O)O)O)O',
      'cellulose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
      'hemicellulose': 'CC(C(C1=CC=CC=C1)O)NC.CC(C(C1=CC=CC=C1)O)NC',
      'pectin': 'C1(C(C(OC(C1O)O)C(=O)O)O)O',
      'xylan': 'C1C(C(C(C(O1)O)O)O)O',
      'arabitol': 'C(C(C(C(CO)O)O)O)O',
      'erythritol': 'C(C(C(CO)O)O)O',
      'galactitol': 'C(C(C(C(C(CO)O)O)O)O)O',
      'glycerol': 'C(C(CO)O)O',
      'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
      'ribitol': 'C(C1C(C(C(O1)C(C(C(C(CO)O)O)O)O)O)O)O',
      'sorbitol': 'C(C(C(C(C(CO)O)O)O)O)O',
      'xylitol': 'C(C(C(C(CO)O)O)O)O',
      'd-minus bornesitol': 'OC1C(O)C(O)C(OC)C(O)C1O',
      'plus inositol': 'C1(C(C(C(C(C1O)O)O)O)O)O',
      'myo inositol': 'C1(C(C(C(C(C1O)O)O)O)O)O',
      'plus quebrachitol': 'COC1C(C(C(C(C1O)O)O)O)O',
      'galactosamine': 'C(C1C(C(C(C(O1)O)N)O)O)O',
      'glucosamine': 'C(C1C(C(C(C(O1)O)N)O)O)O',
    }

    return smiles

  @staticmethod
  def get_smarts():

    smarts = {
    }

    return smarts
