#!/usr/bin/env python3
#
# GlobalChemExtensions - Amino Acid Adapter

# -------------------------------------------

class AminoAcidConverter(object):

    __version__ = '0.0.1'

    def convert_amino_acid_sequence_to_smiles(self, sequence):

        '''

        Convert the Amino Acid Sequence to SMILES

        Arguments:
            Sequence (String): Sequence of the Amino Acid
        Returns:
            Peptide Backbone (String): Smiles string that is not canonical but follows a pattern.

        '''

        amino_acids_sequence = {
            "A": "C",
            "R": "CCCCNC(N)=N",
            "N": "CC(N)=O",
            "D": "CC(O)=O",
            "B": "CC(O)=O",
            "C": "CS",
            "E": "CCC(O)=O",
            "Q": "CCC(N)=O",
            "Z": "CCC(N)=O",
            "G": "[H]",
            "H": "CC1=CNC=N1",
            "I": "C(CC)([H])C",
            "L": "CC(C)C",
            "K": "CCCCN",
            "M": "CCSC",
            "F": "CC1=CC=CC=C1",
            "P": "C2CCCN2",
            "S": "CO",
            "T": "C(C)([H])O",
            "W": "CCC1=CNC2=C1C=CC=C2",
            "Y": "CC1=CC=C(O)C=C1",
            "V": "C(C)C"
        }

        peptide_n_terminus = 'NCC('
        peptide_c_terminus = 'O)=O)'

        for i in range(len(sequence)):

            peptide_n_terminus += 'NCC('
            peptide_c_terminus += '=O)'

        peptide_backbone = peptide_n_terminus + peptide_c_terminus[:-1]

        for i in range(1, len(sequence) + 1):

            start_index_pattern = peptide_backbone.find('NCC')
            peptide_backbone = peptide_backbone[0:start_index_pattern + 1] + 'C([*:' + str(i) + '])' + peptide_backbone[start_index_pattern + 2:]

        for i in range(0, len(sequence)):

            letter = sequence[i]
            smiles_to_add = amino_acids_sequence[letter]

            peptide_backbone = peptide_backbone.replace('[*:%s]' % (i + 1), smiles_to_add)


        return peptide_backbone

    def convert_smiles_to_amino_acid_sequence(self, smiles):

        '''

        Convert Smiles to Amino Acid Sequence

        '''

        amino_acids_sequence = {
            "C" :"A",
            "CCCCNC(N)=N":"R",
            "CC(N)=O":"N",
            "CC(O)=O":"D",
            "CS": "C",
            "CCC(O)=O":"E",
            "CCC(N)=O":"Q",
            "[H]":"G",
            "CC1=CNC=N1" :"H",
            "C(CC)([H])C" :"I",
            "CC(C)C" :"L",
            "CCCCN" :"K",
            "CCSC" :"M",
            "CC1=CC=CC=C1" :"F",
            "C2CCCN2" :"P",
            "CO" :"S",
            "C(C)([H])O" :"T",
            "CCC1=CNC2=C1C=CC=C2" :"W",
            "CC1=CC=C(O)C=C1":"Y",
            "C(C)C" :"V",
        }

        import re

        pattern = re.compile('NC\(.*?\)C\(', flags=re.MULTILINE)
        matches = re.findall(pattern, smiles)

        sequence = ''

        for match in matches:

            match = match[3:-3]

            sequence += amino_acids_sequence[match]

        return sequence