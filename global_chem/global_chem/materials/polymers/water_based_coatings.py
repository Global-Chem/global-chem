#!/usr/bin/env python3
#
# GlobalChem - Water based coatings
#
# ---------------------------------------

class WaterBasedCoatings(object):

    def __init__(self):

        self.name = 'water_based_coatings'

    @staticmethod
    def get_smiles():
        smiles = {
        'n,n-dimethylacrylamide': 'CN(C)C(=O)C=C',
        'acrylamide': 'C=CC(N)=O',
        'acrylic acid,': 'C=CC(=O)O',
        'methacrylic acid': 'CC(=C)C(=O)O',
        'styrene': 'C=Cc1ccccc1',
        'acrylonitrile': 'C=CC#N',
        'methyl methacrylate': 'C=CO',
        'butyl acrylate': 'CCCCOC(=O)C=C',
        'ethyl acrylate': 'CCOC(=O)C=C',
        '2-ethyl hexyl acrylate': 'CCCCC(CC)COC(=O)C=C',
        'hexamethylene diisocyanate': 'C(CCCN=C=O)CCN=C=O',            # Think there is a mistake here
        'toluene diisocyanate': 'CC1=CC=CC=C1.C(=[N-])=O.C(=[N-])=O',  # Think there is a mistake here
        'polyethylene glycol': 'COCCCCO',
        'biuret': 'C(=O)(N)NC(=O)N',
        'vinyl acetate': 'CC(=O)OC=C',
        'veova-10': 'CC(C)(C)CCCCCC(=O)OC=C',
      }

        
        return smiles

