class MangoVitamins(object):

    def __init__(self):

        self.name = 'mango_vitamins'

    @staticmethod
    def get_smiles():

        smiles = {
            'ascorbic acid ': 'C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O',
            'thiamine': 'CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO',
            'riboflavin': 'CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@@H]([C@@H]([C@@H](CO)O)O)O',
            'niacin': 'C1=CC(=CN=C1)C(=O)O',
            'pantothenic acid': 'CC(C)(CO)[C@H](C(=O)NCCC(=O)O)O',
            'pyridoxine ': 'CC1=NC=C(C(=C1O)CO)CO',
            'folic acid ': 'C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N',
            'vitamin A': 'CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C',
            'vitamin E': 'CC1=C(C2=C(CC[C@@](O2)(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=C1O)C)C',
            'vitamin K': 'CC1=C(C(=O)C2=CC=CC=C2C1=O)C/C=C(\C)/CCCC(C)CCCC(C)CCCC(C)C',
        }

        return smiles
