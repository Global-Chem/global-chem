#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem RNA

# -------------------------------------------

# RDKit Imports
# -------------

from rdkit import Chem
from dna_features_viewer import GraphicFeature, GraphicRecord

class GlobalChemRNA(object):

    __RNA__NUCLEOTIDES = {
        'A': 'C1=NC2=NC=NC(=C2N1)N',
        'U': 'C1=CNC(=O)NC1=O',
        'G': 'C1=NC2=C(N1)C(=O)NC(=N2)N',
        'C': 'C1=C(NC(=O)N=C1)N'
    }
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

    def initialize_record(self):

        '''

        Initialize the Graphic Record of the RNA Sequence

        '''

        self.record = GraphicRecord(
            sequence= self.rna_sequence,
            features = self.features
        )

    def label_feature(self, start, end, label='feature', color='#ffcccc'):

        '''

        Features that the user would like to label

        Arguments:
            start (Int): start strand of the rna
            end (Int): end strand of the rna
            label (String): Label of the strand
            color (String): colour of the label

        '''

        feature = GraphicFeature(
            start = start,
            end = end,
            strand = +1,
            color= color,
            label = label
        )

        self.features.append(feature)

    def visualize_strand(self):

        '''

        Visualize the RNA Strand

        Return:
            ax (matplotlib Object): matplotlib Object of the RNA strand.

        '''

        self.initialize_record()

        ax, _ = self.record.plot(figure_width=5)
        self.record.plot_sequence(ax)
        self.record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})

        return ax

    def save_to_image(self):

        '''

        Save the RNA visualization strand to an image

        '''

        self.initialize_record()

        ax, _ = self.record.plot(figure_width=5)
        self.record.plot_sequence(ax)
        self.record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})

        ax.figure.savefig('sequence_and_translation.png', bbox_inches='tight')

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