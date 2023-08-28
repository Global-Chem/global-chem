#!/usr/bin/env python3
#
# GlobalChem - Thai Spices
#
# -----------------------------------------

class ThaiSpices(object):

    def __init__(self):

        self.name = 'thai_spices'

    @staticmethod
    def get_smiles():

        smiles = {
            "3-caren-5-one": "CC1=CC(=O)C2C(C1)C2(C)C",
            
        }

        return smiles

    @staticmethod
    def get_smarts():

        smarts = {
           
        }

        return smarts