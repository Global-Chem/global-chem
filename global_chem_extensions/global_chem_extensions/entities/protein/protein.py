#!/usr/bin/env python3
#
# GlobalChemExtensions - GlobalChem Protein

# -----------------------------------------

# Imports
# -------

import pypdb
from biopandas.pdb import PandasPdb

# RDKit Imports
# -------------

from rdkit import Chem
import rdkit.Chem.Descriptors as Descriptors

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

        # PDB File Handling

        self.pdb_file = pdb_file
        self.fetch_pdb = fetch_pdb
        self.peptide_sequence = peptide_sequence

        # Protein Object Handling

        self.protein = None
        self.protein_sequence = None
        self.protein_dataframe = None

        self.smiles_sequence = None
        self.smarts_sequence = None

        # Ligand Object Handling

        self.ligand_name = None
        self.ligand_smiles = None
        self.ligand_dataframe = None

        if self.pdb_file:
            self.protein, self.protein_sequence = self._read_pdb_file()

        if self.fetch_pdb:
            self.protein, self.protein_sequence = self._fetch_pdb()

        if self.peptide_sequence:
            self.protein_sequence = self.peptide_sequence
            self.convert_to_smiles()

        if self.protein:
            try:
                self.protein_dataframe = self.protein.df['ATOM']
                self.ligand_dataframe = self.protein.df['HETATM']
            except Exception as e:
                print ("WARNING: Ligand Not Found in PDB")
                print (e)
                pass

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

    def determine_bostrom_ligand_criteria(self, verbose=False):

        '''

        Does it Match the Drug Like Filter for a PDB File:

        Arguments:
            verbose (Bool): whether the user would like the output to be verbose

        Returns:
            criteria_met (Bool): a boolean to determine if the ligand passed the criteria or not

        Reference:
            1.) Boström, Jonas, et al. “Do Structurally Similar Ligands Bind in a Similar Fashion?” Journal of Medicinal Chemistry, vol. 49, no. 23, Nov. 2006, pp. 6716–25.

        Filter Criteria:

            1.) 80 < molecular weight (Da) < 750
            2.) 10 < number of nonhydrogen atoms < 70
            3.) Must not contain atoms of types other than
                H, C, O, N, F, P, S, Cl, Br, or I
            4.) must contain at least one
                non-carbon/non-hydrogen atom
            5.) must not contain two or more phosphorus atoms
            6.) must not have more than 10 rotatable bonds
            7.) must not be a nucleic acid
            8.) must not be composed only from
                non-lead-like PDB-HET-groupsb
            9.) must not be covalently bound
            10.) must not have protein contacts from the crystal
                packing environment in less than 3 Å distance
            11.) must have contacts with protein in less than
                7Å distance

        '''

        criteria_met = False

        # Fetch the SMILES

        smiles = self.convert_ligand_to_smiles().strip()
        rdkit_molecule = Chem.MolFromSmiles(smiles)

        # Check 1

        molecular_weight = Descriptors.ExactMolWt(rdkit_molecule)

        if 750 > molecular_weight > 80:

            if verbose:
                print ("Passed Check 1 Molecular Weight: %s " % molecular_weight)

        else:

            if verbose:
                print ("Failed Check 1")

            return criteria_met

        # Check 2

        non_hydrogen_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(rdkit_molecule)

        if 70 > non_hydrogen_atoms > 10:

            if verbose:
                print ("Passed Check 2 Non Hydrogen Atoms: %s " % non_hydrogen_atoms)

        else:

            if verbose:
                print ("Failed Check 2")

            return criteria_met

        # Check 3

        element_boundaries = ['H', 'C', 'O', 'N', 'F', 'P', 'S', 'Cl', 'Br', 'I']

        atom_symbols = []
        atoms = Chem.rdchem.Mol.GetAtoms(rdkit_molecule)

        for atom in atoms:

            symbol = atom.GetSymbol()
            atom_symbols.append(symbol)

            if symbol not in element_boundaries:

                if verbose:
                    print ("Failed Check 3")

                return criteria_met

        if verbose:
            print ("Passed Check 3 All Atoms Are Within Element Boundaries")

        # Check 4

        if 'H' in atom_symbols:

            atom_symbols = list(filter(lambda x: x != 'H', atom_symbols))

        if 'C' in atom_symbols:

            atom_symbols = list(filter(lambda x: x != 'C', atom_symbols))

        if len(atom_symbols) > 0:

            if verbose:
                print ("Passed Check 4 Non-Hydrogen & Non-Carbon Atoms Present: %s" % len(atom_symbols))

        else:

            if verbose:
                print ("Failed Check 4")
            return criteria_met

        # Check 5

        phosphorous_atoms = list(filter(lambda x: x == 'P', atom_symbols))

        if len(phosphorous_atoms) > 1:

            if verbose:
                print ("Failed Check 5")
                return criteria_met

        else:

            if verbose:
                print ("Passed Check 5 Less than Two Phosphorous atoms Present")

        # Check 6

        rotatable_bonds = Descriptors.NumRotatableBonds(rdkit_molecule)

        if rotatable_bonds > 10:

            if verbose:
                print ("Failed Check 6")
                return criteria_met

        else:

            if verbose:
                print ("Passed Check 6 Less than Ten Rotatable Bonds Present")

        # Check 7

        nucleic_acid_smiles = 'NC1OC(COP(O)(O)=O)C(O)C1'
        nucleic_acid_molecule = Chem.MolFromSmiles(nucleic_acid_smiles)

        substructures = rdkit_molecule.GetSubstructMatches(nucleic_acid_molecule)

        if substructures:
            if verbose:
                print ("Failed Check 7")
            return criteria_met

        else:
            if verbose:
                print ("Passed Check 7 No Nucleic Acid Template Found")

        # Check 8

        # if verbose:
        #     print ("Passed Check 8 PDB-HET Groups Not available")

        # Check 9

        if 'COVALENT' in self.protein.pdb_text or 'COVALNT' in self.protein.pdb_text:
            print ("Failed Check 9")
            return criteria_met

        else:
            if verbose:
                print ("Passed Check 9 Not a Covalent Inhibitor")

        # Check 10
        #
        # if verbose:
        #     print ("Passed Check 9 Crystal Packing Contacts Not Implemented")

        # Check 11

        for index, atom_row in self.ligand_dataframe.iterrows():

            # Check to see if there is water

            if atom_row['residue_name'] == 'HOH':
                continue

            x_coord = atom_row['x_coord']
            y_coord = atom_row['y_coord']
            z_coord = atom_row['z_coord']

            reference_point = (x_coord, y_coord, z_coord)

            distances = self.protein.distance(xyz=reference_point, records=('ATOM',))
            atoms_within_seven_angstroms = self.protein.df['ATOM'][distances < 7.0]

            if len(atoms_within_seven_angstroms) > 0:

                print ("Check 11 Has Contacts within 7 Angstroms")
                criteria_met = True
                return criteria_met

        if verbose:
            if not criteria_met:
                print ("Check 11 Failed")

    def convert_ligand_to_smiles(self):


        '''

        Converts the PDB Ligand to SMILES

        '''

        chain_id = self.ligand_dataframe['residue_name'].to_list()[0]

        chem_desc = pypdb.get_info(chain_id,  url_root = 'https://data.rcsb.org/rest/v1/core/chemcomp/')
        chem_desc = pypdb.to_dict(chem_desc)

        self.ligand_name = chem_desc['chem_comp']['name']

        ligand_descriptors = chem_desc['pdbx_chem_comp_descriptor']
        for software in ligand_descriptors:
            if software['program'] == 'OpenEye OEToolkits' and software['type'] == 'SMILES':
                self.ligand_smiles = software['descriptor']

        return self.ligand_smiles