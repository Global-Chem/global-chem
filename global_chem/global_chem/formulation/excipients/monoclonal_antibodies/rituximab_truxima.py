#!/usr/bin/env python3
#
# GlobalChem - Rituximab Truxima
#
# -----------------------------------

class Rituximab_Truxima(object):

    def __init__(self):

        self.name = 'rituximab_truxima'

    @staticmethod
    def get_smiles():

        smiles = {
            'sodium citrate dihydrate': 'C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O.O.O.[Na+].[Na+].[Na+]',
            'polysorbate 80': 'CCCCCCCC/C=C/CCCCCCCC(=O)OCC(C1C(C(CO1)O)O)O',
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
            'sodium citrate dihydrate': '',
            'polysorbate 80': '',
        }

        return smarts
