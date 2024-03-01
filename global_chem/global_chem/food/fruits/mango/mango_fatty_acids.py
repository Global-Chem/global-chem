#!/usr/bin/env python3
#
# GlobalChem - MangoFattyAcids
# -------------------------------

class MangoFattyAcids(object):
    
    '''
    
    Missing Entries
    
                'alpha-Linoleic acid': '',

    '''
    
    def __init__(self):

        self.name = 'mango_fatty_acids'

    @staticmethod
    def get_smiles():

        smiles = {
            'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
            'stearic acid ': 'CCCCCCCCCCCCCCCCCC(=O)O',
            'arachidic acid': 'CCCCCCCCCCCCCCCCCCCC(=O)O',
            'lignoceric acid': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)O',
            'oleic acid': r'CCCCCCCC/C=C\CCCCCCCC(=O)O',
            'linoleic acid': r'CCCCC/C=C\C/C=C\CCCCCCCC(=O)O',
            'myristic acid': 'CCCCCCCCCCCCCC(=O)O',
            'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O ',
            'behenic acid': 'CCCCCCCCCCCCCCCCCCCCCC(=O)O',
            'palmitoleic acid': 'CCCCCC/C=C\CCCCCCCC(=O)O',
            'hexadecenoic acid': 'CCCC/C=C/CCCCCCCCCC(=O)O',
            'heptadecenoic acid': 'CCCCCC/C=C/CCCCCCCCC(=O)O',
            'octadecenoic acid': 'CCCCCC/C=C/CCCCCCCCCC(=O)O',
            'eicosenoic acid': 'CCCCCCCC/C=C\CCCCCCCCCC(=O)O',
            '9,12-Hexadecadienoic acid': 'CCC/C=C/C/C=C/CCCCCCCC(=O)O',
            '9,15-Octadecadienoic acid': 'CC/C=C/CCCC/C=C/CCCCCCCC(=O)O',
            'hepta-2,4(E,E)-dienoic acid': '',
            'linolenic acid': 'CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O'
        }

        return smiles
    
    
    @staticmethod
    def get_smarts():

      smarts = {
              
          }
      
      return smarts
