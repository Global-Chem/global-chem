#!/usr/bin/env python3
#
# GlobalChem - Thermoplastics
#
# ---------------------------------------

class Thermoplastics(object):

    def __init__(self):

        self.name = 'thermoplastics'

    @staticmethod
    def get_smiles():
        smiles = {
        'bisphenol a': 'CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O',
        'bisphenol f': 'C1=CC(=CC=C1CC2=CC=C(C=C2)O)O',
        'chloromethyloxirane': 'C1C(O1)CCl',
        'formaldehyde': 'C=O',
        'epsilon caprolactone': 'C1CCC(=O)OCC1',
        'propylene oxide': 'CC1CO1',
        'ethylene oxide': 'C1CO1',
        'urea': 'C(=O)(N)N',
        'phenol': 'C1=CC=C(C=C1)O',
        'melamine': 'C1(=NC(=NC(=N1)N)N)N',
      }

        
        return smiles

