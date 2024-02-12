#!/usr/bin/env python3
#
# GlobalChem - Solvent Based Coatings
#
# ---------------------------------------

class SolventBasedCoatings(object):

    def __init__(self):

        self.name = 'solvent_based_coatings'

    @staticmethod
    def get_smiles():
        smiles = {
        'phathlic anhydride': 'C1=CC=C2C(=C1)C(=O)OC2=O',
        'isopthalic acid': 'C1=CC(=CC(=C1)C(=O)O)C(=O)O',
        'terepthalic acid': 'C1=CC(=CC=C1C(=O)O)C(=O)O',
        'diethylene glycol': 'C(COCCO)O',
        'adipic acid': 'C(CCC(=O)O)CC(=O)O',
        'pentanedioyl dichloride': 'C(CC(=O)Cl)CC(=O)Cl',
        'benzene-1,4-diol': 'C1=CC(=CC=C1O)O',
        'ricinoleic acid': 'CCCCCCC(CC=CCCCCCCCC(=O)O)O',
        'oleic acid': r'CCCCCCCC/C=C\CCCCCCCC(=O)O',
        'linoleic acid': 'CCCCCC=CCC=CCCCCCCCC(=O)O',
        'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O',
        'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
      }

        
        return smiles

