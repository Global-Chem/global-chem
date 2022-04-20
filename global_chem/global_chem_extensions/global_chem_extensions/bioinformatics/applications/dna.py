#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem DNA

# -------------------------------------------

# RDKit Imports
# -------------

from rdkit import Chem
from dna_features_viewer import GraphicFeature, GraphicRecord

class GlobalChemDNA(object):

    __DNA__NUCLEOTIDES = {
        'A': 'C1=NC2=NC=NC(=C2N1)N',
        'T': 'CC1=CNC(=O)NC1=O',
        'G': 'C1=NC2=C(N1)C(=O)NC(=N2)N',
        'C': 'C1=C(NC(=O)N=C1)N'
    }
    __version__ = '0.0.1'

    def __init__(
            self,
            dna_sequence,
            name = None
    ):

        '''

        Arguments:
            dna_sequence (String): DNA Sequence strand
            name (String): Name of the DNA object if the user wants it.

        '''

        self.name = name
        self.record = None
        self.features = []
        self.dna_sequence = dna_sequence
        self.smiles = ''
        self.sugar_backbone = ''
        self.sequence_length = len(dna_sequence)
        self.dna_molecule = None

    def initialize_record(self):

        '''

        Initialize the Graphic Record of the DNA Sequence

        '''

        self.record = GraphicRecord(
            sequence= self.dna_sequence,
            features = self.features
        )

    def label_feature(self, start, end, label='feature', color='#ffcccc'):

        '''

        Features that the user would like to label

        Arguments:
            start (Int): start strand of the dna
            end (Int): end strand of the dna
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

        Visualize the DNA Strand

        Return:
            ax (matplotlib Object): matplotlib Object of the DNA strand.

        '''

        self.initialize_record()

        ax, _ = self.record.plot(figure_width=5)
        self.record.plot_sequence(ax)
        self.record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})

        return ax

    def save_to_image(self):

        '''

        Save the DNA visualization strand to an image

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

        Convert the DNA sequence to SMILES

        '''

        self.create_sugar_backbone_layer()
        self.smiles = self.sugar_backbone

        for i in range(self.sequence_length):

            modified_molecule = Chem.ReplaceSubstructs(
                Chem.MolFromSmiles(self.smiles),
                Chem.MolFromSmiles('[*:' + str( i + 1 ) + ']'),
                Chem.MolFromSmiles(self.__DNA__NUCLEOTIDES[self.dna_sequence[i]])
            )

            self.smiles = Chem.MolToSmiles(modified_molecule[0], isomericSmiles=False, canonical=True)

        return self.smiles

    def convert_to_smarts(self):

        '''

        Convert the RNA sequence to SMARTS

        '''

        self.convert_to_smiles()

        return Chem.MolToSmarts(Chem.MolFromSmiles(self.smiles))
