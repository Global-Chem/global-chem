#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Mother Nature

# ------------------------------------------------

class BrocaArea(object):

    __version__ = '0.0.1'

    def __init__(self, smiles, iupac, chinese_iupac):

        '''

        Arguments:
            smiles (String): smiles string to translate
            iupac (String): iupac string name to translate
            chinese_iupac (String): chinese iupac name

        '''

        self.smiles = smiles
        self.iupac = iupac
        self.chinese_iupac = chinese_iupac

    def translate_smiles_to_iupac(self):

        '''

        Translate the SMILES to IUPAC

        '''

        pass




