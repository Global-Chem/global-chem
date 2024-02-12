#!/usr/bin/env python3
#
# GlobalChem - Rubbers
#
# ---------------------------------------

class Rubbers(object):

    def __init__(self):

        self.name = 'rubber'

    @staticmethod
    def get_smiles():
       
        smiles = {
          '2-chlorobuta-1,3-diene': 'C=CC(=C)Cl',
          'isoprene': 'CC(=C)C=C',
          'ethylidene norbornene': 'CC=C1CC2CC1C=C2',
          'dicyclopentadiene': 'C1C=CC2C1C3CC2C=C3',
          'vinyl norbornene': 'C=CC1CC2CC1C=C2',
          'vinylidene fluoride': 'C=C(F)F',
          'ethylene': 'C=C',
          'butadiene': 'C=CC=C',
        }

        return smiles