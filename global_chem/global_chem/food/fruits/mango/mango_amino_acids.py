class MangoAminoAcids(object):

    def __init__(self):

        self.name = 'mango_amino_acids'

    @staticmethod
    def get_smiles():

        smiles = {
            'arginine': 'C(CC(C(=O)O)N)CN=C(N)N',
            'cysteine': 'C([C@@H](C(=O)O)N)S',
            'glycine': 'C(C(=O)O)N',
            'isoleucine': 'CC[C@H](C)[C@@H](C(=O)O)N ',
            'leucine': 'CC(C)C[C@@H](C(=O)O)N',
            'lysine': 'C(CCN)C[C@@H](C(=O)O)N',
            'methionine': 'CSCC[C@@H](C(=O)O)N',
            'phenylalanine': 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N',
            'proline': 'C1C[C@H](NC1)C(=O)O',
            'serine': 'C([C@@H](C(=O)O)N)O',
            'threonine': 'C[C@H]([C@@H](C(=O)O)N)O',
            'tryptophan': 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N ',
            'tyrosine': 'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O  ',
            'valine': 'CC(C)[C@@H](C(=O)O)N',
            'glutamic acid': 'C(CC(=O)O)[C@@H](C(=O)O)N',
            'aspartic acid': 'C([C@@H](C(=O)O)N)C(=O)O',
            'alanine': 'C[C@@H](C(=O)O)N',
            'arginine': 'C(C[C@@H](C(=O)O)N)CN=C(N)N',
            'histidine': 'C1=C(NC=N1)C[C@@H](C(=O)O)N'
        }

        return smiles


    @staticmethod
    def get_smarts():

      smarts = {
              
          }
      
      return smarts