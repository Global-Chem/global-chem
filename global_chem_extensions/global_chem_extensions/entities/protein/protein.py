#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Protein

# -----------------------------------------

# Imports
# -------

from rdkit import Chem
from biopandas.pdb import PandasPdb

class GlobalChemProtein(object):


    __version__ = '0.0.1'

    def __init__(
            self,
            pdb_file = None,
            fetch_pdb = None,
            peptide_sequence = None,
    ):
        '''

        Arguments:
            pdb_file (String): path to pdb file
            fetch_pdb (String): PDB ID
            peptide_sequence (String): Reading the protein peptide sequence

        '''

        self.pdb_file = pdb_file
        self.fetch_pdb = fetch_pdb
        self.peptide_sequence = peptide_sequence

        self.protein = None
        self.protein_sequence = None
        self.smiles_sequence = None
        self.smarts_sequence = None

        if self.pdb_file:
            self.protein, self.protein_sequence = self._read_pdb_file()

        if self.fetch_pdb:
            self.protein, self.protein_sequence = self._fetch_pdb()

        if self.peptide_sequence:
            self.protein_sequence = self.peptide_sequence
            self.convert_to_smiles()

    def _read_pdb_file(self):

        '''

        Use Biopandas to Read the PDB

        Returns:
            biopandas_pdb (Biopandas PDB Object): The biopandas object
            sequence (String): the peptide sequence string

        '''

        biopandas_pdb = PandasPdb().read_pdb(self.pdb_file)

        sequence_dataframe = biopandas_pdb.amino3to1()
        sequence = ''.join(sequence_dataframe['residue_name'].to_list())

        return biopandas_pdb, sequence

    def _fetch_pdb(self):

        '''

        Use Biopandas to Fetch the PDB File

        Returns:
            biopandas_pdb (Biopandas PDB Object): The biopandas object
            sequence (String): the peptide sequence string

        '''

        biopandas_pdb = PandasPdb().fetch_pdb(self.fetch_pdb)

        sequence_dataframe = biopandas_pdb.amino3to1()
        sequence = ''.join(sequence_dataframe['residue_name'].to_list())

        return biopandas_pdb, sequence

    def convert_to_smiles(
            self,
            mark_nitrogen_backbone = False,
            mark_carbonyl_carbon_backbone = False,
            # mark_oxygen = False,
    ):

        '''

        Convert to SMILES

        Arguments:
            mark_nitrogen_backbone (Bool): whether to mark the nitrogen backbone
            mark_carbonyl_carbon_backbone (Bool): whether to mark the carbonyl carbon on the backbone
            mark_oxygen (Bool): whether to mark the double bond oxygen


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

        n_terminus_pattern = 'NCC('
        c_terminus_pattern = 'O)=O)'

        replacement_counter = 0

        if mark_nitrogen_backbone:
            replacement_counter += 1
            n_terminus_pattern = n_terminus_pattern.replace('N', '[N:1]')

        if mark_carbonyl_carbon_backbone:
            replacement_counter += 1
            n_terminus_pattern = n_terminus_pattern.replace('C(', '[C:%s](' % replacement_counter)

        # if mark_oxygen:
        #     replacement_counter += 1
        #     c_terminus_pattern = c_terminus_pattern.replace('=O', '(=[O:%s])' % replacement_counter)

        peptide_n_terminus = n_terminus_pattern
        peptide_c_terminus = c_terminus_pattern

        for i in range(len(self.protein_sequence)):

            peptide_n_terminus += n_terminus_pattern

            # if mark_oxygen:
            #     peptide_c_terminus = '([O:%s])' % replacement_counter
            # else:
            peptide_c_terminus += '=O)'

        peptide_backbone = peptide_n_terminus + peptide_c_terminus[:-1]

        for i in range(1, len(self.protein_sequence) + 1):

            start_index_pattern = peptide_backbone.find(n_terminus_pattern)

            if mark_nitrogen_backbone:
                peptide_backbone = peptide_backbone[0:start_index_pattern + 5] + 'C([*:' + str(i) + '])' + peptide_backbone[start_index_pattern + 6:]
            else:
                peptide_backbone = peptide_backbone[0:start_index_pattern + 1] + 'C([*:' + str(i) + '])' + peptide_backbone[start_index_pattern + 2:]

        for i in range(0, len(self.protein_sequence)):

            letter = self.protein_sequence[i]
            smiles_to_add = amino_acids_sequence[letter]

            peptide_backbone = peptide_backbone.replace('[*:%s]' % (i + 1), smiles_to_add)

        self.smiles_sequence = peptide_backbone

        return self.smiles_sequence

    def convert_to_smarts(
            self,
            mark_nitrogen_backbone = False,
            mark_carbonyl_carbon_backbone = False,
            # mark_oxygen = False,
    ):

        """

        Convert to a Smarts String

        """

        smiles = self.convert_to_smiles(
            mark_nitrogen_backbone = mark_nitrogen_backbone,
            mark_carbonyl_carbon_backbone = mark_carbonyl_carbon_backbone,
        )

        self.smarts_sequence = Chem.MolToSmarts(Chem.MolFromSmiles(smiles))

        return self.smarts_sequence

if __name__ == '__main__':

    gc_protein = GlobalChemProtein(
        # pdb_file='file.pdb',
        # fetch_pdb='5tc0',
        peptide_sequence='AAAA',
    )

    print (gc_protein.convert_to_smarts(
        mark_nitrogen_backbone=False,
        mark_carbonyl_carbon_backbone = False,
    ))
