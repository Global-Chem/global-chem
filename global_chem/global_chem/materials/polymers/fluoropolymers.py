#!/usr/bin/env python3
#
# GlobalChem - FluoroPolymers
#
# ---------------------------------------

class FluoroPolymers(object):

    def __init__(self):

        self.name = 'fluoropolymers'

    @staticmethod
    def get_smiles():
        smiles = {
          'perfluorocycloalkene': 'C1(=C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)F)C(F)(F)F',
          'vinyl fluoride': 'C=CF',
          'vinylidene': 'C=C(F)F',
          'tetrafluoroethylene': 'C(=C(F)F)(F)F',
          'chlorotrifluoroethylene': 'C(=C(F)Cl)(F)F',
          'hexafluoropropylene': 'C(=C(F)F)(C(F)(F)F)F',
        }

        
        return smiles

