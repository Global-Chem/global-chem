class MangoFattyAcids(object):

    def __init__(self):

        self.name = 'mango_fatty_acids'

    @staticmethod
    def get_smiles():

        smiles = {
            'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
            'stearic acid ': 'CCCCCCCCCCCCCCCCCC(=O)O',
            'arachidic acid': 'CCCCCCCCCCCCCCCCCCCC(=O)O',
            'lignoceric acid': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)O',
            'oleic acid': 'CCCCCCCC/C=C\CCCCCCCC(=O)O',
            'linoleic acid': 'CCCCC/C=C\C/C=C\CCCCCCCC(=O)O',
            'alpha-Linoleic acid': '',
            'myristic acid': 'CCCCCCCCCCCCCC(=O)O',
            'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O',
            'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O ',
            'arachidic acid': 'CCCCCCCCCCCCCCCCCCCC(=O)O',
            'behenic acid': 'CCCCCCCCCCCCCCCCCCCCCC(=O)O',
            'lignoceric acid': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)O',
            'palmitoleic acid': 'CCCCCC/C=C\CCCCCCCC(=O)O',
            'hexadecenoic acid': 'CCCC/C=C/CCCCCCCCCC(=O)O',
            'heptadecenoic acid': 'CCCCCC/C=C/CCCCCCCCC(=O)O',
            'oleic acid': 'CCCCCCCC/C=C\CCCCCCCC(=O)O',
            'octadecenoic acid': 'CCCCCC/C=C/CCCCCCCCCC(=O)O',
            'eicosenoic acid': 'CCCCCCCC/C=C\CCCCCCCCCC(=O)O',
            '9,12-Hexadecadienoic acid': 'CCC/C=C/C/C=C/CCCCCCCC(=O)O',
            'linoleic acid': 'CCCCC/C=C\C/C=C\CCCCCCCC(=O)O',
            '9,15-Octadecadienoic acid': 'CC/C=C/CCCC/C=C/CCCCCCCC(=O)O',
            'hepta-2,4(E,E)-dienoic acid': '',
            'linolenic acid': 'CC/C=C\C/C=C\C/C=C\CCCCCCCC(=O)O'
        }

        return smiles
