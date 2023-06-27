class MangoPhenolicAcids(object):

    def __init__(self):

        self.name = 'mango_phenolic_acids'

    @staticmethod
    def get_smiles():

        smiles = {
            'gallic Acid': 'C1=C(C=C(C(=C1O)O)O)C(=O)O',
            'vanillic acid ': 'COC1=C(C=CC(=C1)C(=O)O)O',
            'syringic acid ': 'COC1=CC(=CC(=C1O)OC)C(=O)O',
            'protocatechuic acid': 'C1=CC(=C(C=C1C(=O)O)O)O',
            'para hydroxybenzoic acid ': 'C1=CC(=CC=C1C(=O)O)O',
            'paracoumaric acid': 'C1=CC(=CC=C1/C=C/C(=O)O)O',
            'chlorogenic acid ': 'C1[C@H]([C@H]([C@@H](C[C@@]1(C(=O)O)O)OC(=O)/C=C/C2=CC(=C(C=C2)O)O)O)O',
            'ferulic acid': 'COC1=C(C=CC(=C1)/C=C/C(=O)O)O',
            'caffeic acid': 'C1=CC(=C(C=C1/C=C/C(=O)O)O)O',
            'theogallin': 'C1[C@H]([C@H]([C@@H](C[C@@]1(C(=O)O)O)OC(=O)C2=CC(=C(C(=C2)O)O)O)O)O',
        }

        return smiles
