#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Amylopectin

# ---------------------------------------------

# RDKit Imports
# -------------

from rdkit import Chem

class GlobalChemAmylopectin(object):

    __version__ = '0.0.1'

    def __init__(
            self,
            rna_sequence,
            name = None
    ):

        '''

        Arguments:
            rna_sequence (String): RNA Sequence strand
            name (String): Name of the RNA object if the user wants it.

        '''

        self.name = name
        self.record = None
        self.features = []
        self.rna_sequence = rna_sequence
        self.smiles = ''
        self.sugar_backbone = ''
        self.sequence_length = len(rna_sequence)
        self.rna_molecule = None

    def create_sugar_backbone_layer(self):

        '''

        Creates the sugar backbone layer

        '''

        counter = 2

        start = '[*:1]C1CC(OP(=O)([O-])[O-])C'

        if self.sequence_length == 1:

            self.sugar_backbone = start + '(CO)O1'
            return self.smiles

        start = '[*:1]C1CC(OP(=O)([O-])[O-])C(COP(=O)([O-])'

        self.sugar_backbone += start

        for j in range(self.sequence_length - 1):

            self.sugar_backbone += 'OC2CC([*:%s])OC2COP(=O)([O-])' % counter
            counter += 1

        end = 'OC2CCOC2CO)O1'

        self.sugar_backbone += end

    def convert_to_smiles(self):

        '''

        Convert the RNA sequence to SMILES

        '''

        self.create_sugar_backbone_layer()
        self.smiles = self.sugar_backbone

        for i in range(self.sequence_length):

            modified_molecule = Chem.ReplaceSubstructs(
                Chem.MolFromSmiles(self.smiles),
                Chem.MolFromSmiles('[*:' + str( i + 1 ) + ']'),
                Chem.MolFromSmiles(self.__RNA__NUCLEOTIDES[self.rna_sequence[i]])
            )

            self.smiles = Chem.MolToSmiles(modified_molecule[0], isomericSmiles=False, canonical=True)

        return self.smiles

    def convert_to_smarts(self):

        '''

        Convert the RNA sequence to SMARTS

        '''

        self.convert_to_smiles()

        return Chem.MolToSmarts(Chem.MolFromSmiles(self.smiles))