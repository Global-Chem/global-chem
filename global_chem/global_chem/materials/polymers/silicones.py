#!/usr/bin/env python3
#
# GlobalChem - Silicones
#
# ---------------------------------------

class Silicones(object):

    def __init__(self):

        self.name = 'silicones'

    @staticmethod
    def get_smiles():
        smiles = {

          'octamethylcyclotetrasiloxane': 'C[Si]1(O[Si](O[Si](O[Si](O1)(C)C)(C)C)(C)C)',
        }

        
        return smiles

